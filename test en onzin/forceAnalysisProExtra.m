
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Top level section 2: FX2_LOAD_RAW  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1
    
    %%% This is the second section in a series of routines to analyze
    %%% force-extension data taken with the new multi bead code for MT
    %%% (J. P. Cnossen, D. Dulin and N. H. Dekker Rev Sci Instrum 85, (2014).)
    %%%
    %%% This section is used to load the data for the force-extension
    %%% measurements, namely the large matrix of position traces for different
    %%% magnet heights.
    %%%
    %%% The script assumes that you have determined the z-offsets previously and
    %%% loads them z-offsets from a file
    %%%
    %%% Author: Jan Lipfert
    %%% Date: 2013-09-16
    
    clear all; clc; close all;
    
    
    %%% Index for each data set
    DATA = 1;
    
    
    %%% Write here some descriptive text about your data
    %%% ---
    %%% Example of force extension data
    %%% 21 kbp DNA, M270 beads, 1 mm gap verticaly oriented magnets
    if DATA ==1
%         traces_file = 'bead.txt';
%         motors_file = 'bead_motors.txt';
%        zoffsets_file = 'FX_offsets.txt';
        freq = 144; %%% Acquisition frequency in Hz
    end
    
    
    
    
    
    %%%--- Read in the previously determined z-offsets from file
    zoffsets = 0; %nm
    loaddata = load('matlab.mat');
    names = fieldnames(loaddata);
    
    
    %%%--- Read in data ---
    dataBead = double(loaddata.(names{1}));
    data = dataBead(:,1:3);
    zmag = loaddata.RdRp4(:,4);
    
    Nbeads = 1;
    
    %%%--- Parse the bead data ---
    for i=1:Nbeads
        display(['Parsing data for bead ' num2str(i) ' of ' num2str(Nbeads)])
        bead(i).time = 1:length(data(:,1));
        bead(i).x = data(:,2)*1000; %nm
        bead(i).y = data(:,1)*1000;  %%%x data is y data, want x en y zijn omgewisseld in daldrop
        bead(i).z = data(:,3)*1000;
        
    end
    
    
    %%% Subtract previously determined z-offsets
    if 1
        if length(zoffsets) == Nbeads
            for i=1:Nbeads
                bead(i).z = bead(i).z - zoffsets(i);
            end
            display('Subtracted pre-determined z-offsets.')
        else
            display(['Number of beads in the big data file (' num2str(Nbeads) ') does not agree with the number of beads in the zoffset file (' num2str(length(zoffsets)) ')'])
        end
    end
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Top level section 3:  FX3_ANALYZE %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if 1
    
    %%% This is the third section in a series of routines to analyze
    %%% force-extension data taken with the new multi bead code for MT
    %%% (J. P. Cnossen, D. Dulin and N. H. Dekker Rev Sci Instrum 85, (2014).)
    %%%
    %%% This section assumes that you have read the (z-offset corrected) raw
    %%% data traces
    %%%
    %%% Author: Jan Lipfert
    %%% Date: 2013-09-16, updated 2014-07-17
    
    
    clc;
    
    plotflag = 0; %%% Whether or not to show the fits
    
    
    %%% Some variables ----------------------------------------------
    kT = 4.1; %pN nm
    eta = 10E-10; %viscosity in pN s/nm^2
    Rbead = 1400; % Bead radius in nm
    nBlock = 5;
    
    %%% ---------------------------------------------------------------
    %%% Figure out where the magnets are moving, from the motor file
    %%% This determines the "plateaus", where the nmagnet height is
    %%% constant and where we we want to analyze the forces
    %%% ---------------------------------------------------------------
    Nsmooth_zmag = 200;
    Nsmooth_dzmag = 200;
    small = 10^(-5); %%% Threshold to determine where it is moving
    
    %%% Some variables for error checking
    Nmin_points_plat = 100; % Minimum number of points in a plateau
    zmags = [];
    
    if 1
        zmag_smooth = smooth(zmag, Nsmooth_zmag, 'moving');
        diff_zmag = smooth(diff(zmag_smooth), Nsmooth_dzmag, 'lowess');
        
        %%% Find all the "plateaus", i.e. the sets of points where the magnets
        %%% are not moving up or down
        Nfirst =1; Nlast  =1;
        
        %%% Check whether we are starting in a plateau
        if (abs(diff_zmag(1)) < small & abs(diff_zmag(2)) < small)
            tplat(1).first = 1;
            Nfirst = Nfirst +1;
        end
        
        for i=2:length(diff_zmag)
            if (abs(diff_zmag(i)) < small & abs(diff_zmag(i-1)) > small)
                tplat(Nfirst).first = i;
                Nfirst = Nfirst +1;
            end
            
            if (abs(diff_zmag(i)) > small & abs(diff_zmag(i-1)) < small)
                tplat(Nlast).last = i;
                Nlast = Nlast +1;
            end
        end
        
        %%% Check whether we are ending in a plateau
        if (abs(diff_zmag(end)) < small & abs(diff_zmag(end-1)) < small)
            tplat(Nlast).last = length(diff_zmag);
            Nlast = Nlast +1;
        end
        
        Nplat = Nfirst-1;
        display(['Found ' num2str(Nplat) ' Zmag plateaus']);
        %display(['Nfirst ' num2str(Nfirst) ' Nlast ' num2str(Nlast)]);
        
        
        %%% Throw out plateaus that have fewer than Nmin_points_plat points
        %%% This is not quite it yet, since it keeps the plateaus separate that
        %%% are interupted by a "fake" (i.e. too short) plateau, but getting it
        %%% perfect is tricky
        if 1
            count = 1;
            Ngoodplat = 0;
            for j = 1:Nplat
                if length(tplat(j).first:tplat(j).last) > Nmin_points_plat
                    plat(count).first = tplat(j).first;
                    plat(count).last = tplat(j).last;
                    
                    Ngoodplat = Ngoodplat + 1;
                    count = count + 1;
                    
                else
                    display(['Plateau # ' num2str(j) ' has only ' num2str(length(Nmin_points_plat)) ' points!' ])
                end
            end
            
            Nplat = Ngoodplat;
        end
        
        
        
        %%% Plot the magnet height information, after smoothing, with
        %%% plateaus annotated
        if 1
            figure(1);clf; hold on; box on; %%% Magnet rot. vs. time
            platind = find(abs(diff_zmag) < small);
            faketime = 1:length(zmag_smooth);
            plot(faketime , zmag_smooth, 'b-')
            plot(faketime(platind), zmag_smooth(platind), 'r.')
            for j = 1:Nplat
                plot(faketime(plat(j).first) , zmag_smooth(plat(j).first), 'ko', 'markersize', 10)
                plot(faketime(plat(j).last) , zmag_smooth(plat(j).last), 'mo', 'markersize', 10)
                %%% Capture the output
                zmags = [zmags zmag_smooth(plat(j).last)];
            end
            set(gca, 'fontsize', 14, 'fontweight', 'bold', 'linewidth', 1,'TickLength',[0.02 0.02])
            xlabel('Time'); ylabel('Zmag')
            title(['Red = plateaus; Black / magneta circles = start / stop; Found ' num2str(Nplat) ' plateaus.' ])
            
            figure(2);clf; hold on; box on; %%% Derivative of magnet rot. vs. time
            plot(1:length(diff_zmag) , diff_zmag, 'b-')
            set(gca, 'fontsize', 14, 'fontweight', 'bold', 'linewidth', 1,'TickLength',[0.02 0.02])
            xlabel('Time'); ylabel('dzmag/dt')
        end
        
    end
    %%% ------- END PLATEAU FINDIND --------- %%%
    
    pause;
    
    
    
    %%% ------------------------------------------------------- %%%
    %%% Look over the beads and determine the force for each plateau
    %%% ------------------------------------------------------- %%%
    if 1
        for i=1:Nbeads
            
            bead(i).ext = zeros(1,Nplat);
            bead(i).Fx_real = zeros(1,Nplat);
            bead(i).Fy_real = zeros(1,Nplat);
            bead(i).PSDforce = zeros(1,Nplat);
            bead(i).PSD2force = zeros(1,Nplat);
            bead(i).Rfitx = zeros(1,Nplat);
            bead(i).PSD3force = zeros(1,Nplat);
            bead(i).Rfity = zeros(1,Nplat);
            
            
            %%% Use the script "analyze_one_trace" to determine the force for each trace
            for k=1:Nplat
                
                [Ext, Fx_real, Fy_real, PSDfit, PSDforce, MLfitx, MLforcex, Rfitx, MLfity, MLforcey, Rfity]=...
                analyze_one_trace2(bead(i).time(plat(k).first:plat(k).last),...
                    bead(i).x(plat(k).first:plat(k).last),... 
                    bead(i).y(plat(k).first:plat(k).last),...
                    bead(i).z(plat(k).first:plat(k).last),...
                freq, Rbead, kT, eta, nBlock);

                display(['Finished working on bead # ' num2str(i) ', plateau number ' num2str(k)])
                bead(i).ext(k) = Ext;
                bead(i).Fx_real(k) = Fx_real;
                bead(i).Fy_real(k) = Fy_real;
                PSDfit = PSDfit;
                bead(i).PSDforce(k) = PSDforce;
                PSD2fitx = MLfitx;
                bead(i).PSD2force(k) = MLforcex;
                bead(i).Rfitx(k) = Rfitx;
                PSD3fity = MLfity;
                bead(i).PSD3force(k) = MLforcey;
                bead(i).Rfity(k) = Rfity;
                bead(i).difference(k) = abs(PSDforce - MLforcey);

                %pause;
            end
        end
    end
end