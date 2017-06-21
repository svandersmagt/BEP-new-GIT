
% % % Matlab package for calibrating forces in Magnetic Tweezers (MT) measurements
% % % ------------------------------------------------------------------------------------
% % % Toivo Hensgens, Jan Lipfert, and Zhongbo Yu, Nynke Dekker lab, updated on 2014-Dec-3
% % % ====================================================================================
% % %
% % % We used this script to calculate forces in magnetic tweezers. This
% % % script is the same as the one on the website of Nynke Dekker lab:
% % % http://nynkedekkerlab.tudelft.nl/software/.
% % %
% % % --------------
% % % User’s actions
% % % --------------
% % %
% % % 1. You make a new script in Matlab using the following codes. The Matlab version we used is R2013b. 
% % %
% % % 2. You can find the sample data at the website of Nynke Dekker's lab:
% % % http://nynkedekkerlab.tudelft.nl/software/.
% % %  
% % % 3. You can directly run the script with the sample data in the same folder.
% % %
% % % 4. OPTIONAL: Changing 'if 1' to 'if 0' at the sections of the top level, you can block the running of the script to see the intermediate results.
% % %
% % % ----------------------------------------------
% % % Introdcution of this script
% % % ----------------------------------------------
% % %
% % % There are four sections at the top level of this script:
% % % The first section, FX1_DETERMINE_ZOFFSET, loads data from the data source and determine the z-offset, from a measurement at very low force. Here you need find the sample data at http://nynkedekkerlab.tudelft.nl/software/, Or you use your own data. 
% % % The second section, FX2_LOAD_RAW, loads data from the data source.
% Again, you need to prepare the data as mentioned above.
% % % The third section, FX3_ANALYZE, analyzes the forces.
% % % The fourth section, FX4_ANALYZE_Lp_Lc, summarizes the worm-like chain
% fitting results over all the beads analysed. Here, we have just one bead
% in the sample data.
% % %
% % %
% % % Functions at the low level:
% % % The function of analyze_one_trace takes a time trace and calculates the force by calling the 3 different force analysis scripts (see below)
% % %
% % % All the m files other than the top and sub-functions of m files:
% % % Most of these files are self-explanatory. Some important ones are:
% % % The m file of analyze_PSD.m computes the force of a time trace by the PSD method
% % % The m file of analyze_AV.m computes the force of a time trace by the AV method
% % % The m file of spectral_analysis.m computes the force of a time trace by the SA method
% % %
% % %
% % % ----------------
% % % Sample data
% % % ----------------
% % % You need the following smple data to test this script, which you can find at http://nynkedekkerlab.tudelft.nl/software/:
% % % data.txt
% % % data_motors.txt
% % % offset.txt
% % % offset_motors.txt
% % % This script calculates the forces using the coordinates of a sphere in the field of magnets. The coordinates of x, y, and z are generally evenly spaced time traces. In this particular data of the example, we use a custom-writing Labview software (ref 1) to record the traces which is saved in the file of ‘bead.txt’. We limit the size of files by only presenting the bead fluctuations of y and z dimensions, which have been corrected by a reference bead. The file of ‘bead_motors.txt’ contains the z positions of the magnets at the time corresponding to the coordinates of the sphere in ‘bead.txt’. The files of 'offset.txt' and 'offset_motors.txt' contains the bead fluctuations without magnets.
% % %
% % %
% % % ------------------------------
% % % Explanation of jargons
% % % ------------------------------
% % % F_IND:  a variable to select the dimension of fluctuation for calculating forces, i.e., x trace (1) or y trace (2). In the example, we present the trace of y only.
% % % AV: Allan variance, see ref 2
% % % SA: Spectral analysis, see ref 2
% % % PSD: power spectra density, see ref 3
% % % Ext: DNA extension
% % % WLC: Worm-like chain model, see ref 4
% % %
% % %
% % %
% % %
% % % -------------------------
% % % Acknowledgement
% % % -------------------------
% % % We have used the m file allan v3.0 by Fabian Czerwinski, exp2fit by Per Sundqvist, from the file exchange in matlab central
% % %
% % %
% % % ---------------
% % % References:
% % % ---------------
% % % 1.Cnossen, J. P., et al. (2014). "An optimized software framework for real-time, high-throughput tracking of spherical beads." Rev Sci Instrum 85(10): 103712.
% % % 2.Lansdorp BM & Saleh OA (2012) Power spectrum and Allan variance methods for calibrating single-molecule video-tracking instruments. Rev Sci Instrum 83(2):025115.
% % % 3.te Velthuis AJ, Kerssemakers JW, Lipfert J, & Dekker NH (2010) Quantitative guidelines for force calibration through spectral analysis of magnetic tweezers data. Biophys J 99(4):1292-1302.
% % % 4.Bouchiat, et al. Biophys J 76:409 (1999)



function forceAnalysis



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Top level section 1: FX1_DETERMINE_ZOFFSET  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    
    %%% This is the first section in a series of routines to analyze
    %%% force-extension data taken with the new multi bead code for MT
    %%% (J. P. Cnossen, D. Dulin and N. H. Dekker Rev Sci Instrum 85, (2014).)
    %%%
    %%% This section is used to determine the z-offset for each bead
    %%% ---
    %%% To determine the z-offsets, record bead position traces under conditions
    %%% where the beads are essentially at the surface, either by going to very
    %%% low (essentially zero) stretching force or by going to a large number
    %%% of supercoils /turns under low forces.
    %%%
    %%% - Set the file names and path informaiton in the "DATA" block below.
    %%% The variable "DATA" controls which data sets are read and analyzed.
    %%%
    %%% - The code determines the offset by smoothing the position traces with a
    %%% window of Nsmooth points abd takes the minimum of the smoothed traces
    %%% as the surface position = z-offset
    %%%
    %%% - The results get saved into a file of which the name is set through the
    %%% "output_name" variable. Format is first colum = bead number, second
    %%% column = offset in mum
    %%%
    %%% Author: Jan Lipfert
    %%% Date: 2013-09-16
    
    clear all; clc; close all;
    
    
    %%% Index to select which data set to run
    DATA = 1;
    
    
    %%% Write here some descriptive text about your data
    %%% ---
    %%% Example data set on FX measurements; M270 beads, 21 kbp DNA
    if DATA ==1
        traces_file = '..\BEPmaarnietgit\offset.txt';
        motors_file = '..\BEPmaarnietgit\offset_motors.txt';
        output_name = '..\BEPmaarnietgit\FX_offsets3.txt';
        Nref = 1;
    end
    
    Nsmooth = 100;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%--- Read in data ---
    %%%%%%%%%%%%%%%%%%%%%%%
    
    data = load(traces_file);
    zmag = load(motors_file); % for this example, motor file contains just the data of zmag
    
    
    %%%--- Parse the bead data ---
    % For this example, there is just one bead, which is already corrected by reference
    % bead. The dimensions of y and z are the only columns saved for this example.
    for i=1
        display(['Parsing data for bead ' num2str(i) ' of ' num2str(1)])
        bead(i).time = 1:length(data(:,1));
        bead(i).y = data(:,1);
        bead(i).z = data(:,2);
    end
    
    %%% --- Plot the motor information --- %%%
    if 1
        %%% Plot the magnet height information
        figure(100);clf; hold on; box on; %%% Magnet rot. vs. time
        faketime = 1:length(zmag);
        plot(faketime , zmag, 'r-')
        set(gca, 'fontsize', 16, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02]);
        xlabel('Time (frames)'); ylabel('Zmag (mm)')
        title('Motor height')
        
    end
    
    z_offsets = [];
    %%% Loop over beads and show the height information
    if 1
        
        
        for i=1 % for this example, there is just one bead
            figure(1); clf; hold on; box on;
            plot(bead(i).time, bead(i).z, 'k-', 'linewidth', 1)
            
            %%% Smooth and find minimum
            smooth_z = smooth(bead(i).z, Nsmooth, 'moving');
            plot(bead(i).time, smooth_z, 'r-', 'linewidth', 2)
            
            [min_z ind] = min(smooth_z);
            plot(bead(i).time(ind), smooth_z(ind), 'bx', 'linewidth', 2, 'markersize', 20)
            
            z_offsets = [z_offsets min_z];
            
            set(gca, 'fontsize', 16, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02]);
            xlabel('Time (s)'); ylabel('z (um)')
            title(['Bead # ' num2str(i)])
            pause;
            
        end
    end
    
    
    %%% Plot the offsets
    if 1
        figure(2); clf; hold on; box on;
        plot(1, z_offsets, 'bo', 'linewidth', 2, 'markersize', 20)
        set(gca, 'fontsize', 16, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02]);
        xlabel('Bead #'); ylabel('z-offset (um)')
        
    end
    
    %%% Save the data
    if 1
        foo = [(1)' z_offsets'];
        save(output_name, 'foo', '-ascii')
        
    end
    
    
end



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
        traces_file = '..\BEPmaarnietgit\fx_11_19Nov2013.txt';
        motors_file = '..\BEPmaarnietgit\fx_11_19Nov2013_motors.txt';
        zoffsets_file = '..\BEPmaarnietgit\FX_offsets3.txt';
        Nref = 1;
        F_IND = 2; %%% whether to use x or y (1=x, 2=y)
        freq = 60; %%% Acquisition frequency in Hz
    end
    
    
    
    
    
    %%%--- Read in the previously determined z-offsets from file
    zoff_data = load(zoffsets_file);
    zoffsets = zoff_data(:,2);
    
    
    %%%--- Read in data ---
    data = load(traces_file);
    zmag = load(motors_file);
    zmag = zmag(:,3);
    
    Nbeads = 1;
    
    %%%--- Parse the bead data ---
    for i=1:Nbeads
        display(['Parsing data for bead ' num2str(i) ' of ' num2str(Nbeads)])
        bead(i).time = 1:length(data(:,1));
        bead(i).y = data(:,1);
        bead(i).z = data(:,2);
        
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
    
    plotflag = 1; %%% Whether or not to show the fits
    
    
    %%% Some variables ----------------------------------------------
    kT = 4.1;
    Rbead = 1.4;      % Bead radius in mum
    
    
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
            bead(i).AVforce = zeros(1,Nplat);
            bead(i).fcorner = zeros(1,Nplat);
            bead(i).SAforce = zeros(1,Nplat);
            bead(i).fcorner2 = zeros(1,Nplat);
            
            
            %%% Use the script "analyze_one_trace" to determine the force for each trace
            for k=1:Nplat
                if length(plat(k).first:plat(k).last) > Nmin_points_plat
                    
                    [bead(i).ext(k),...
                        bead(i).Fx_real(k),...
                        bead(i).Fy_real(k),...
                        bead(i).PSDforce(k),...
                        PSDfit,...
                        bead(i).AVforce(k),...
                        AVfit,...
                        bead(i).fcorner(k),...
                        SAfit,...
                        bead(i).SAforce(k),...
                        bead(i).fcorner2(k)] = ...
                        analyze_one_trace(bead(i).time(plat(k).first:plat(k).last),...
                        bead(i).y(plat(k).first:plat(k).last),...
                        bead(i).y(plat(k).first:plat(k).last),...
                        bead(i).z(plat(k).first:plat(k).last),...
                        i, F_IND, freq, plotflag); %For this example, we replace x dimension with y dimension to keep the size of data small
                    
                    display(['Working on bead # ' num2str(i) ', plateau number ' num2str(k)])
                    
                    if plotflag == 1
                        %                     pause;
                    end
                end
            end
            
            
            %%% Fit a two exponential model to the data
            %%% The model: F0 + Amp * exp(Zmag/l_decay) + Amp * exp(Zmag/l_decay)
            %%% Note that this is just an empirical fit, no real physics
            s=exp2fit(zmags, bead(i).PSDforce,2);
            bead(i).PSD_F0 = s(1);
            bead(i).PSD_Amp = s(2);
            bead(i).PSD_2_decay = s(3);
            bead(i).PSD_Amp_2 = s(4);
            bead(i).PSD_2_decay_2 = s(5);
            bead(i).PSD_Zfit = min(zmags):0.01:max(zmags);
            bead(i).PSD_Ffit = s(1)+s(2)*exp(-bead(i).PSD_Zfit/s(3))+s(4)*exp(-bead(i).PSD_Zfit/s(5));
            
            s=exp2fit(zmags, bead(i).AVforce,2);
            bead(i).AV_F0 = s(1);
            bead(i).AV_Amp = s(2);
            bead(i).AV_2_decay = s(3);
            bead(i).AV_Amp_2 = s(4);
            bead(i).AV_2_decay_2 = s(5);
            bead(i).AV_Zfit = min(zmags):0.01:max(zmags);
            bead(i).AV_Ffit = s(1)+s(2)*exp(-bead(i).AV_Zfit/s(3))+s(4)*exp(-bead(i).AV_Zfit/s(5));
            
            s=exp2fit(zmags, bead(i).SAforce,2);
            bead(i).SA_F0 = s(1);
            bead(i).SA_Amp = s(2);
            bead(i).SA_2_decay = s(3);
            bead(i).SA_Amp_2 = s(4);
            bead(i).SA_2_decay_2 = s(5);
            bead(i).SA_Zfit = min(zmags):0.01:max(zmags);
            bead(i).SA_Ffit = s(1)+s(2)*exp(-bead(i).SA_Zfit/s(3))+s(4)*exp(-bead(i).SA_Zfit/s(5));
            
            if 1
                %%% Fit the WLC chain model to the data
                %%% This is real physics and determines Lp (persistence length)
                %%% and Lc (the contour length)
                [bead(i).PSD_Lp, bead(i).PSD_Lc, bead(i).PSD_Forcefit, bead(i).PSD_Extfit] = inExtWLCfit(bead(i).ext, bead(i).PSDforce, 0.1);
                [bead(i).AV_Lp, bead(i).AV_Lc, bead(i).AV_Forcefit, bead(i).AV_Extfit] = inExtWLCfit(bead(i).ext, bead(i).AVforce, 0.1);
                [bead(i).SA_Lp, bead(i).SA_Lc, bead(i).SA_Forcefit, bead(i).SA_Extfit] = inExtWLCfit(bead(i).ext, bead(i).SAforce, 0.1);
            end
            
            
            %%% Show a Force vs. Zmag plot for this bead
            if 1
                
                %%% Plot F vs. Zmag, linear
                figure(3); clf; hold on; box on;
                plot(zmags, bead(i).Fx_real, 'ko', 'markersize', 15, 'linewidth', 1)
                plot(zmags, bead(i).PSDforce, 'bo', 'markersize', 15, 'linewidth', 2)
                plot(zmags, bead(i).AVforce, 'ro', 'markersize', 15, 'linewidth', 2)
                plot(zmags, bead(i).SAforce, 'go', 'markersize', 15, 'linewidth', 2)
                
                %%% Plot the fits
                plot(bead(i).PSD_Zfit,  bead(i).PSD_Ffit, 'b-', 'linewidth', 2)
                plot(bead(i).AV_Zfit,  bead(i).AV_Ffit, 'r-', 'linewidth', 2)
                plot(bead(i).SA_Zfit,  bead(i).SA_Ffit, 'g-', 'linewidth', 2)
                
                
                set(gca, 'fontsize', 18, 'fontweight', 'bold', 'linewidth', 1,'TickLength',[0.02 0.02])
                xlabel('Magnet position (mm)');
                ylabel('Force (pN)');
                legend('Real space', 'PSD', 'Allan', 'SA')
                %%% Show the fitted parameter from the Allan variance fit in the title
                title(['From AV fit: Amp=' num2str(bead(i).AV_Amp,3) ';  l(decay)=' num2str(bead(i).AV_2_decay,3) '; Amp2=' num2str(bead(i).AV_Amp_2,3) ';  l(decay)2=' num2str(bead(i).AV_2_decay_2,3) ';   F0=' num2str(bead(i).AV_F0,2) ])
                
                %%% Plot F vs. Zmag, log scaling
                figure(4); clf; hold on; box on;
                plot(zmags, bead(i).Fx_real, 'ko', 'markersize', 15, 'linewidth', 1)
                plot(zmags, bead(i).PSDforce, 'bo', 'markersize', 15, 'linewidth', 2)
                plot(zmags, bead(i).AVforce, 'ro', 'markersize', 15, 'linewidth', 2)
                plot(zmags, bead(i).SAforce, 'go', 'markersize', 15, 'linewidth', 2)
                
                %%% Plot the fits
                plot(bead(i).PSD_Zfit,  bead(i).PSD_Ffit, 'b-', 'linewidth', 2)
                plot(bead(i).AV_Zfit,  bead(i).AV_Ffit, 'r-', 'linewidth', 2)
                plot(bead(i).SA_Zfit,  bead(i).SA_Ffit, 'g-', 'linewidth', 2)
                
                set(gca, 'fontsize', 18, 'fontweight', 'bold', 'linewidth', 1,'TickLength',[0.02 0.02])
                xlabel('Magnet position (mm)');
                ylabel('Force (pN)');
                set(gca, 'yscale', 'log')
                legend('Real space', 'PSD', 'Allan', 'SA')
                title(['Bead # ' num2str(i)])
                
            end
            
            
            
            
            %%% Show a FX plot for this bead
            if 1
                
                %%% Plot force vs. DNA extension
                figure(5); clf; hold on; box on;
                plot(bead(i).ext, bead(i).Fx_real, 'ko', 'markersize', 15, 'linewidth', 1)
                
                plot(bead(i).ext, bead(i).PSDforce, 'bo', 'markersize', 15, 'linewidth', 2)
                plot(bead(i).ext, bead(i).AVforce, 'ro', 'markersize', 15, 'linewidth', 2)
                plot(bead(i).ext, bead(i).SAforce, 'go', 'markersize', 15, 'linewidth', 2)
                
                if 1
                    plot(bead(i).PSD_Extfit, bead(i).PSD_Forcefit, 'b-', 'markersize', 15, 'linewidth', 2)
                    plot(bead(i).AV_Extfit, bead(i).AV_Forcefit, 'r-', 'markersize', 15, 'linewidth', 2)
                    plot(bead(i).SA_Extfit, bead(i).SA_Forcefit, 'g-', 'markersize', 15, 'linewidth', 2)
                end
                
                set(gca, 'fontsize', 18, 'fontweight', 'bold', 'linewidth', 1,'TickLength',[0.02 0.02])
                xlabel('Extension (um)');
                ylabel('Force (pN)');
                legend('Real space', 'PSD', 'Allan', 'SA', 'location', 'northwest')
                title(['From AV fit: Lp=' num2str(bead(i).AV_Lp*1000,3) 'nm;  Lc=' num2str(bead(i).AV_Lc,3) 'um'  ])
                
                
                %%% Log plot
                figure(6); clf; hold on; box on;
                plot(bead(i).ext, bead(i).Fx_real, 'ko', 'markersize', 15, 'linewidth', 1)
                
                plot(bead(i).ext, bead(i).PSDforce, 'bo', 'markersize', 15, 'linewidth', 2)
                plot(bead(i).ext, bead(i).AVforce, 'ro', 'markersize', 15, 'linewidth', 2)
                plot(bead(i).ext, bead(i).SAforce, 'go', 'markersize', 15, 'linewidth', 2)
                
                if 1
                    plot(bead(i).PSD_Extfit, bead(i).PSD_Forcefit, 'b-', 'markersize', 15, 'linewidth', 2)
                    plot(bead(i).AV_Extfit, bead(i).AV_Forcefit, 'r-', 'markersize', 15, 'linewidth', 2)
                    plot(bead(i).SA_Extfit, bead(i).SA_Forcefit, 'g-', 'markersize', 15, 'linewidth', 2)
                end
                
                set(gca, 'fontsize', 18, 'fontweight', 'bold', 'linewidth', 1,'TickLength',[0.02 0.02])
                xlabel('Extension (um)');
                ylabel('Force (pN)');
                
                set(gca, 'yscale', 'log')
                legend('Real space', 'PSD', 'Allan', 'SA', 'location', 'northwest')
                title(['Bead # ' num2str(i)])
            end
            
            
            %%% Show a corner freq vs. Zmag plot for this bead
            if 1
                figure(7); clf; hold on; box on;
                plot(zmags, bead(i).fcorner2, 'ro', 'markersize', 15, 'linewidth', 2)
                
                set(gca, 'fontsize', 14, 'fontweight', 'bold', 'linewidth', 1,'TickLength',[0.02 0.02])
                xlabel('Magnet position (mm)');
                ylabel('f(corner, Hz)');
                
            end
            
            %         pause;
            
        end
    end
    
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Top level section 4: FX4_ANALYZE_Lp_Lc %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    
    %%% This is the forth section in a series of routines to analyze
    %%% force-extension data taken with the new multi bead code for MT
    %%% (J. P. Cnossen, D. Dulin and N. H. Dekker Rev Sci Instrum 85, (2014).)
    %%%
    %%% This section assumes that you have carried out fits of the WLC model for
    %%% all beads and all magnet heights/forces
    %%%
    %%% Author: Jan Lipfert
    %%% Date:  2014-07-17
    
    
    
    
    %%% Collect the results of the WLC fits, for the various force
    %%% determination methods
    PSD_Lps = [];
    PSD_Lcs = [];
    AV_Lps = [];
    AV_Lcs = [];
    SA_Lps = [];
    SA_Lcs = [];
    
    for i=1:Nbeads
        
        PSD_Lps = [PSD_Lps bead(i).PSD_Lp];
        PSD_Lcs = [PSD_Lcs bead(i).PSD_Lc];
        AV_Lps = [AV_Lps bead(i).AV_Lp];
        AV_Lcs = [AV_Lcs bead(i).AV_Lc];
        SA_Lps = [SA_Lps bead(i).SA_Lp];
        SA_Lcs = [SA_Lcs bead(i).SA_Lc];
        
    end
    
    
    %%% Histogram the results
    if 1
        centersVals = [0:5:200];
        
        figure(8); clf; box on; hold on;
        subplot(1,2,1); box on; hold on;
        [counts centers] = hist(PSD_Lps*1000, centersVals);
        plot(centers, counts, 'bo-','linewidth', 2)
        
        [counts centers] = hist(AV_Lps*1000, centersVals);
        plot(centers, counts, 'ro-','linewidth', 2)
        
        [counts centers] = hist(SA_Lps*1000, centersVals);
        plot(centers, counts, 'go-','linewidth', 2)
        
        set(gca,'LineWidth', 1,'FontSize', 15, 'FontWeight', 'bold', 'TickLength',[0.02 0.02])
        xlabel('Lp (nm)'); ylabel('Counts')
        title(['Lp(PSD) = ' num2str(median(PSD_Lps*1000),3) ' nm'  ])
        
        
        centersVals = [0:0.2:8];
        
        subplot(1,2,2); box on; hold on;
        [counts centers] = hist(PSD_Lcs, centersVals);
        plot(centers, counts, 'bo-','linewidth', 2)
        
        [counts centers] = hist(AV_Lcs, centersVals);
        plot(centers, counts, 'ro-','linewidth', 2)
        
        [counts centers] = hist(SA_Lcs, centersVals);
        plot(centers, counts, 'go-','linewidth', 2)
        
        set(gca,'LineWidth', 1,'FontSize', 15, 'FontWeight', 'bold', 'TickLength',[0.02 0.02])
        xlabel('Lc (um)'); ylabel('Counts')
        title(['Lp(PSD) = ' num2str(median(PSD_Lcs),3) ' um' ])
        legend('PSD', 'Allan', 'Aartjan', 'location', 'northwest')
        
    end
    
    
end
end





























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Low level function  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Ext, Fx_real, Fy_real, PSDforce, PSDfit, AVforce, AVfit, fcorner, SAfit, SAforce, fcorner2] = analyze_one_trace(time, x, y, z, beadnumber, F_IND, fs, plotflag)
%
% Function to analyze magnetic tweezers time traces
%
% Input: (time, x, y, z, filenr, F_IND, fs, plotflag)
% - time in s (vector)
% - trace x in mum (vector, same length)
% - trace y in mum (vector, same length)
% - trace z in mum (vector, same length)
% - beadnumber (integer; display in the title of the (x,y,z) plot)
% - F_IND - set 1 to analyze x fluctuations, set 2 to analyze y fluctuations
% - fs in Hz (number, sampling frequency)
% if plotflag == 1: Show the x,y,z time trace and histograms, as well as
% the fits for the PSD, AV and SA methods

foo = 1;

%%% Global, "hard-coded" variables
binWidth_XY = 0.01; % in mum
binWidth_Z  = 0.01; % in mum
kT = 4/1000;

%%% Analyze the real space fluctuations
mean_x = mean(x);
mean_y = mean(y);
mean_z = mean(z);
std_x  = std(x);
std_y  = std(y);
std_z  = std(z);

Fx_real = kT*mean_z./std_x^2;
Fy_real = kT*mean_z./std_y^2;
Ext = mean_z;

%%% Analyze the fluctuations in freq domain using Power Spectral Density (PSD) max likelihood fit
%%% Also analyze the real space fluctuations using the Allan variance
%%% These implement the two methods of Landsdorp and Saleh (RSI 2012)
%%% The third method is that by te Velthuis et al. (Bioph J 2010)

if F_IND == 1;
    [PSDfit, PSDforce, fcorner] = analyze_PSD(fs,mean_z,x,plotflag);
    [AVfit, AVforce] = analyze_AV(fs,mean_z,x,plotflag);
    [SAfit, SAforce, fcorner2] = spectral_analysis(fs,mean_z,x,plotflag);
elseif F_IND == 2;
    [PSDfit, PSDforce, fcorner] = analyze_PSD(fs,mean_z,y,plotflag);
    [AVfit, AVforce] = analyze_AV(fs,mean_z,y,plotflag);
    [SAfit, SAforce, fcorner2] = spectral_analysis(fs,mean_z,y,plotflag);
end


if plotflag == 1
    %%% If plotflag = 1, show the x,y,z-traces with histograms
    %%% and the x-y, x-z and y-z plots
    %%% Code from show_xyz
    warning off
    figure(5);clf;
    
    subplot(3,3,1); hold on; box on;
    title([ 'F_x = ' num2str(Fx_real, 2) ' pN'], 'fontsize', 14);
    plot(time, x)
    ylabel('X (\mum)', 'fontsize', 14);
    xlabel('Time  (s)', 'fontsize', 14);
    set(gca,'LineWidth', 1,'FontSize', 14)
    set(gca,'TickLength',[0.02 0.02])
    
    subplot(3,3,2); hold on; box on;
    title([ 'F_y = ' num2str(Fy_real, 2) ' pN'], 'fontsize', 14);
    plot(time, y, 'r')
    ylabel('Y (\mum)', 'fontsize', 14);
    xlabel('Time  (s)', 'fontsize', 14);
    set(gca,'LineWidth', 1,'FontSize', 14)
    set(gca,'TickLength',[0.02 0.02])
    
    
    subplot(3,3,3); hold on; box on;
    title([ '<z> = ' num2str(mean_z, 3) ' \mum'], 'fontsize', 14);
    plot(time, z, 'Color', [0 0.5 0])
    ylabel('Z (\mum)', 'fontsize', 14);
    xlabel('Time  (s)', 'fontsize', 14);
    set(gca,'LineWidth', 1,'FontSize', 14)
    set(gca,'TickLength',[0.02 0.02])
    
    
    %%% Plot X vs. Y
    subplot(3,3,7); hold on; box on;
    title([ 'Bead number: ' num2str(beadnumber)], 'fontsize', 14);
    
    plot(x, y, 'k')
    
    maxval = max([max(x) abs(min(x)) max(y) abs(min(y))]);
    
    if maxval > 2.1
        axis([-3 3 -3 3])
    elseif maxval > 1.1
        axis([-2 2 -2 2])
    else
        axis([-1 1 -1 1])
    end
    axis square
    xlabel('X (\mum)', 'fontsize', 14);
    ylabel('Y (\mum)', 'fontsize', 14);
    set(gca,'LineWidth', 1,'FontSize', 14)
    set(gca,'TickLength',[0.02 0.02])
    
    %%% Plot X vs. Z
    subplot(3,3,8); hold on; box on;
    plot(x, z, 'k')
    
    xlabel('X (\mum)', 'fontsize', 14);
    ylabel('Z (\mum)', 'fontsize', 14);
    set(gca,'LineWidth', 1,'FontSize', 14)
    set(gca,'TickLength',[0.02 0.02])
    
    %%% Plot Y vs. Z
    subplot(3,3,9); hold on; box on;
    plot(y, z, 'k')
    
    xlabel('Y (\mum)', 'fontsize', 14);
    ylabel('Z (\mum)', 'fontsize', 14);
    set(gca,'LineWidth', 1,'FontSize', 14)
    set(gca,'TickLength',[0.02 0.02])
    
    
    
    %%% Subtract out the means
    x = x - mean_x;
    y = y - mean_y;
    z = z - mean_z;
    binWidth = 0.001;
    
    
    %%% Analyze X %%%
    %%%--- Fit an exponential distribution ---%%%
    binCenters = min(x):binWidth:max(x);
    [MUHAT,SIGMAHAT,MUCI,SIGMACI] = normfit(x);
    xfit = min(x):(binWidth*0.1):max(x);
    yfit = length(x)* binWidth * normpdf(xfit, MUHAT, SIGMAHAT);
    
    
    subplot(3,3,4); hold on; box on;        %%% X histogramm
    hist(x, binCenters)
    h1 = findobj(gca,'Type','patch');
    set(h1,'FaceColor', [0 0 0],'EdgeColor','k')
    
    plot(xfit, yfit, 'color', [0 0 1], 'linewidth', 2)
    ylabel('Counts', 'fontsize', 14);
    xlabel('X (\mum)', 'fontsize', 14);
    set(gca,'LineWidth', 1,'FontSize', 14)
    set(gca,'TickLength',[0.02 0.02])
    
    title(['\sigma  = ' num2str(SIGMAHAT*1000,3) ' nm ' ])
    
    %%% Analyze Y %%%
    %%%--- Fit an exponential distribution ---%%%
    binCenters = min(y):binWidth:max(y);
    [MUHAT,SIGMAHAT,MUCI,SIGMACI] = normfit(y);
    xfit = min(y):(binWidth*0.1):max(y);
    yfit = length(y)* binWidth * normpdf(xfit, MUHAT, SIGMAHAT);
    
    
    subplot(3,3,5); hold on; box on;        %%% X histogramm
    hist(y, binCenters)
    h1 = findobj(gca,'Type','patch');
    set(h1,'FaceColor', [0 0 0],'EdgeColor','k')
    
    plot(xfit, yfit, 'color', [1 0 0], 'linewidth', 2)
    ylabel('Counts', 'fontsize', 14);
    xlabel('Y (\mum)', 'fontsize', 14);
    set(gca,'LineWidth', 1,'FontSize', 14)
    set(gca,'TickLength',[0.02 0.02])
    
    title(['\sigma  = ' num2str(SIGMAHAT*1000,3) ' nm ' ])
    
    
    %%% Analyze Z %%%
    %%%--- Fit an exponential distribution ---%%%
    binCenters = min(z):binWidth:max(z);
    [MUHAT,SIGMAHAT,MUCI,SIGMACI] = normfit(z);
    xfit = min(z):(binWidth*0.1):max(z);
    yfit = length(z)* binWidth * normpdf(xfit, MUHAT, SIGMAHAT);
    
    
    subplot(3,3,6); hold on; box on;        %%% X histogramm
    hist(z, binCenters)
    h1 = findobj(gca,'Type','patch');
    set(h1,'FaceColor', [0 0 0],'EdgeColor','k')
    
    plot(xfit, yfit, 'color', [0 0.5 0], 'linewidth', 2)
    ylabel('Counts', 'fontsize', 14);
    xlabel('Z (\mum)', 'fontsize', 14);
    set(gca,'LineWidth', 1,'FontSize', 14)
    set(gca,'TickLength',[0.02 0.02])
    
    title(['\sigma  = ' num2str(SIGMAHAT*1000,3) ' nm ' ])
    %print(gcf,'-r300','-dpng','testtrace.png');
    warning on
end

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Low level function  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [avar]=allan(data, tau)

% Compute various Allan deviations for a constant-rate time series
% [AVAR]=allan(DATA, TAU)
%
% INPUTS:
% DATA should be a struct and has the following fields:
%  DATA.freq    the time series measurements in arb. units
%  DATA.rate    constant rate of time series in (Hz)
%               (Differently from previous versions of allan.m,
%               it is not possible to compute variances for time-
%               stamp data anymore.)
% TAU is an array of the tau values for computing Allan deviations
%
% OUTPUTS:
% AVAR is a struct and has the following fields (for values of tau):
%  AVAR.sig     = standard deviation
%  AVAR.sig2    = Allan deviation
%  AVAR.sig2err = standard error of Allan deviation
%  AVAR.osig    = Allan deviation with overlapping estimate
%  AVAR.osigerr = standard error of overlapping Allan deviation
%  AVAR.msig    = modified Allan deviation
%  AVAR.msigerr = standard error of modified Allan deviation
%  AVAR.tsig    = timed Allan deviation
%  AVAR.tsigerr = standard error of timed Allan deviation
%  AVAR.tau1    = measurement interval in (s)
%  AVAR.tauerr  = errors in tau that might occur because of initial
%  rounding
%
% NOTES:
% Calculations of modified and timed Allan deviations for very long time
% series become very slow. It is advisable to uncomment .msig* and .tsig*
% only after calculations of .sig*, .sig2* and .osig* have been proven
% sufficiently fast.
%
% No pre-processing of the data is performed.
% For constant-rate time series, the deviations are only calculated for tau
% values greater than the minimum time between samples and less than half
% the total time.
%
% versionstr = 'allan v3.0';
% FCz OCT2009
% v3.0  faster and very plain code, no plotting; various Allan deviations
%       can be calculated; script and sample data are availabie on
%       www.nbi.dk/~czerwin/files/allan.zip
%		(Normal, overlapping and modified Allan deviations are calculated in one function,
%		 in strong contrast to MAHs approach of splitting up among various functions. This might be beneficial for individual cases though.)
%
% MAH 2009
% v2.0 and others
%
% FCz OCT2008
% v1.71 'lookfor' gives now useful comments; script and sample data are
%       availabie on www.nbi.dk/~czerwin/files/allan.zip
% v1.7  Improve program performance by mainly predefining matrices outside
%       of loops (avoiding memory allocation within loops); no changes to
%       manual
%
% early program core by Alaa MAKDISSI 2003
% (documentation might be found http://www.alamath.com/)
% revision and modification by Fabian CZERWINSKI 2009
%
% For more information, see:
% [1] Fabian Czerwinski, Andrew C. Richardson, and Lene B. Oddershede,
% "Quantifying Noise in Optical Tweezers by Allan Variance,"
% Opt. Express 17, 13255-13269 (2009)
% http://dx.doi.org/10.1364/OE.17.013255


n=length(data.freq);
jj=length(tau);
m=floor(tau*data.rate);

avar.sig     = zeros(1, jj);
avar.sigerr  = zeros(1, jj);
avar.sig2    = zeros(1, jj);
avar.sig2err = zeros(1, jj);
avar.osig    = zeros(1, jj);
avar.osigerr = zeros(1, jj);
% avar.msig    = zeros(1, jj);
% avar.msigerr = zeros(1, jj);
% avar.tsig    = zeros(1, jj);
% avar.msigerr = zeros(1, jj);

%tic;

for j=1:jj
    % fprintf('.');
    
    D=zeros(1,n-m(j)+1);
    D(1)=sum(data.freq(1:m(j)))/m(j);
    for i=2:n-m(j)+1
        D(i)=D(i-1)+(data.freq(i+m(j)-1)-data.freq(i-1))/m(j);
    end
    
    %standard deviation
    avar.sig(j)=std(D(1:m(j):n-m(j)+1));
    avar.sigerr(j)=avar.sig(j)/sqrt(n/m(j));
    
    %normal Allan deviation
    avar.sig2(j)=sqrt(0.5*mean((diff(D(1:m(j):n-m(j)+1)).^2)));
    avar.sig2err(j)=avar.sig2(j)/sqrt(n/m(j));
    
    %overlapping Allan deviation
    z1=D(m(j)+1:n+1-m(j));
    z2=D(1:n+1-2*m(j));
    u=sum((z1-z2).^2);
    avar.osig(j)=sqrt(u/(n+1-2*m(j))/2);
    avar.osigerr(j)=avar.osig(j)/sqrt(n-m(j));
    
    %     %modified Allan deviation
    %     u=zeros(1,n+2-3*m(j));
    %     z1=D(1:m(j));
    %     z2=D(1+m(j):2*m(j));
    %     for L=1:n+1-3*m(j)
    %         u(L)=(sum(z2-z1))^2;
    %         z1=z1-y(L)+y(L+m(j));
    %         z2=z2-y(L+m(j))+y(L+2*m(j));
    %     end
    %     avar.msigerr(j)=avar.msig(j)/sqrt(n-m(j));
    %     uu=mean(u);
    %     avar.msig(j)=sqrt(uu/2)/m(j);
    %
    %     %timed Allan deviation
    %     avar.tsig(j)=tau(j)*avar.msig(j)/sqrt(3);
    %     avar.tsigerr(j)=avar.tsig(j)/sqrt(n-m(j));
    
    % toc
    
end;

avar.tau1=m/data.rate;
avar.tauerr=tau-avar.tau1;

%toc;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Low level function  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AVmodel = analytical_AV_overdamped_bead(alpha,kappa,taus)
% Calculates the theoretical Allan variance for Brownian motion in
% a harmonic well (Lansdorp and Saleh, RSI 2012)
% Input: alpha in pN s/nm (number)
%        kappa in pN/nm (number)
%        taus in s (vector)
% Output: AVmodel in nm^2 (vector)

kT = 4; %pN nm
fac = alpha./(kappa*taus);

AVmodel = 1 + 2*fac.*exp(-1./fac) - fac/2.*exp(-2./fac) - 3*fac/2;
AVmodel = 2 * kT * fac / kappa .* AVmodel;

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Low level function  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PSDmodel = analytical_PSD_overdamped_bead(alpha,kappa,fs,f)
% Calculates the theoretical Power Spectral Density for Brownian motion
% in a harmonic well (Lansdorp and Saleh, RSI 2012)
% Input: alpha in pN s/nm (number)
%        kappa in pN/nm (number)
%        fs in s (number)
%        f in s (vector)
% Output: PSDmodel in nm^2/Hz (vector)

kT = 4; %pN nm
prefac = 2*kT*alpha/kappa^2;

num = 2*alpha/kappa*fs*(sin(pi*f/fs).^2)*sinh(kappa/(alpha*fs));
denom = cos(2*pi*f/fs)-cosh(kappa/(alpha*fs));

PSDmodel = prefac*(1+num./denom);

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Low level function  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [AVfit AVforce] = analyze_AV(fs,mean_z,x,plotflag)
% Analyzes the Allan Variance of a single trace
% and runs fits....
%
% Input:
% - fs: measurement frequency in Hz
% - mean_z: mean extension in um
% - x: data trace in um
% - plotflag: to plot set 1 - to not plot set 0
%
% Output:
% - AVfit: (1) alpha in pN s/nm (2) kappa in pN/nm (3) good/bad (1/0)
% - AVforce: force according to ML PSD fit in pN

%%% Get an estimate of the force
kT = 4/1000; %pN um
F_est = kT*mean_z/std(x)^2; %pN

%%% Add a line to x to make the length 2^integer, change to nm
x(end+1) = x(1);
x = x*1000; %nm
N = length(x);
logN2 = log(N)/log(2);

%%% Calculate the AV for 100 points, logspace sampled
%%% between ts and (N/2)*ts (so from m=1 point to m=N/2)
ts = 1/fs; %time intervals in s
taus = logspace(log(ts)/log(10),log(N/2*ts)/log(10),100);
allanin.freq = x;
allanin.rate = fs;
allanout = allan(allanin,taus);
AV = allanout.sig2.^2;

%%% Now calc the AV for the octave-sampled points
tausOS = 2.^[0:logN2-1]*ts; %these are the octave-sampled values for tau
allanoutOS = allan(allanin,tausOS);
AV_OS = allanoutOS.osig.^2; %osig is the overlapping AV output of the 'allan' function


%%% Try to fit the octave-sampled Allan Variance using the model
%%% by Lansdorp and Saleh (RSI 2012) using Max Likelihood Fit
alpha_0 = 1E-5; kappa_0 = 4E-4;
[par res ef] = fminsearch(@(par) costfunction_AVfit(N*ts, tausOS, par(1), par(2), AV_OS), [alpha_0 kappa_0]);
alpha = par(1); kappa = par(2); AVfit = [alpha kappa];

AVforce = kappa*1000*mean_z; % in pN


%%% Plot if required
if plotflag == 1
    figure(2); clf; hold on; box on;
    plot(taus, AV, 'b-', 'linewidth', 2)
    plot(taus, allanout.osig.^2, '--', 'linewidth', 2,'color',[0.5 0.5 0])
    plot(tausOS, AV_OS, 'rx', 'linewidth', 2)
    
    %%% Calc and plot the model
    AVmodel = analytical_AV_overdamped_bead(alpha,kappa,taus);
    plot(taus, AVmodel, 'k--', 'linewidth', 2)
    
    legend('AV','Overlapping AV','Octave-sampled points','ML fit','location','southwest')
    legend(gca,'boxoff')
    set(gca, 'fontsize', 14, 'linewidth', 1, 'fontweight', 'bold', 'TickLength',[0.02 0.02]);
    set(gca,'yscale','log')
    set(gca,'xscale','log')
    set(gca,'xlim',[min(taus)/2 max(taus)*2])
    xlabel('t (s)')
    ylabel('Allan variance (nm^2)')
    title(['F_{var} = ' num2str(F_est,2) ' pN,   F_{AV} = ' num2str(AVforce,2) ' pN'])
    %print(gcf,'-r300','-dpng','testAV.png');
end


end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Low level function  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [MLfit MLforce fcorner] = analyze_PSD(fs,mean_z,x,plotflag)
% Analyzes the Power Spectral Density of a single trace
% and runs fits
%
% Input:
% - fs: measurement frequency in Hz
% - mean_z: mean extension in um
% - x: data trace in um
% - plotflag: to plot set 1 - to not plot set 0
%
% Output:
% - MLfit: (1) alpha in pN s/nm (2) kappa in pN/nm (3) good/bad (1/0)
% - MLforce: force according to ML PSD fit in pN
% - fcorner: corner frequency (according to Lorentzian fit) in Hz

%%% Get an estimate of the force and the corner frequency
kT = 4/1000; %pN um
F_est = kT*mean_z/std(x-mean(x))^2; %pN
f_c = calc_fcorner(F_est,mean_z); %Hz
fitgood = f_c < fs/2;

%%% Add a line to x to make the length 2^integer, change to nm
%%% Subtract the mean (offset)
x(end+1) = x(1);
x = x*1000; %nm
%x = x - mean(x);
N = length(x);
logN2 = log(N)/log(2);

%%% Calc PSD, find data points below 1/20 of f_c (the corner frequency)
%%% Also throw away first point, which is merely mean(x)
%%% NOTE: 1/20 is hard-coded, seems to work fine for all data
[f PSD ~] = calc_powersp(x,fs);
f(1) = []; PSD(1) = [];
goodinds = f > f_c(1)/20;

%%% Max likelihood fit to PSD (using goodinds)
%%% According to the model by Lansdorp and Saleh (RSI 2012)
alpha_0 = 1E-5; kappa_0 = 4E-4;
[par res ef] = fminsearch(@(par) costfunction_PSDfit(f(goodinds), par(1), par(2),fs, PSD(goodinds)), [alpha_0 kappa_0]);
alpha = par(1); kappa = abs(par(2)); MLfit = [alpha kappa fitgood];
PSDmodel = analytical_PSD_overdamped_bead(alpha,kappa,fs,f);
MLforce = kappa*1000*mean_z;

%%% Fit a Lorentzian to find the corner frequency (using the 'goodinds')
fcorner_0 =1; Amp_0 = 10^(-2);
[par res ef] = fminsearch(@(par) norm(PSD(goodinds) - (par(1)*((1+(f(goodinds)/par(2)).^(2)).^(-1)))), [Amp_0 fcorner_0]);
fcorner = par(2); Amp  = par(1);
fcorner = abs(fcorner);
Lorentzianfit =  (Amp*(1+(f/fcorner).^(2)).^(-1));
[i ii] = min(abs(f-fcorner));

%%% Plot if required
if plotflag == 1
    figure(1); clf; hold on; box on;
    plot(f,PSD,'b','linewidth',1)
    plot(f, PSDmodel, 'k--', 'linewidth', 2)
    %plot(f, Lorentzianfit, 'g--', 'linewidth', 2)
    plot(fcorner,Lorentzianfit(ii),'x','markersize',14,'color',[0 0.8 0],'linewidth',5)
    %line([f_c f_c],[min(PSD) max(PSD)],'color','r','linewidth',2,'linestyle','--')
    plot(f(~goodinds),PSD(~goodinds),'r','linewidth',1.5)
    
    legend('PSD','ML fit','Theoretical F_{corner}','Non-fitted data','location','southwest')
    legend(gca,'boxoff')
    
    set(gca, 'fontsize', 14, 'linewidth', 1, 'fontweight', 'bold', 'TickLength',[0.02 0.02]);
    xlabel('Frequency (Hz)')
    ylabel('PSD (nm^2/Hz)')
    set(gca,'yscale','log')
    set(gca,'xscale','log')
    set(gca,'xlim',[f(1)/2 f(end)*2])
    title(['F_{var} = ' num2str(F_est,2) ' pN,   F_{PSD} = ' num2str(MLforce,2) ' pN'])
    
    if 0
        %%% Calculate the 'blocked' PSD, divide the data into blocks
        %%% of m points. Try to find a decent automated value of m
        %%% NOTE: mfac is thus hard-coded
        mfac = 3/4;
        m = 2^(round(mfac*logN2));
        b = N/m-1; %number of bins
        for i = 1:b
            xbin = x((i-1)*m/2+1:m+(i-1)*m/2);
            %%% 'Window' the data using a Hann window
            if 1
                Hannwindow = hann(m);
                xbin = xbin.*Hannwindow*sqrt(8/3);
            end
            [fbin PSDbin(:,i) c] = calc_powersp(xbin,fs);
        end
        %%% Average to find PSD_blocked, throw away first point = mean(xbin)
        PSD_blocked = mean(PSDbin');
        fbin(1) = []; PSD_blocked(1) = [];
        plot(fbin, PSD_blocked, '--', 'linewidth', 2,'color',[0.5 0 0.5])
        plot(f, PSDmodel, 'k--', 'linewidth', 2)
        plot(fcorner,Lorentzianfit(ii),'x','markersize',14,'color',[0 0.8 0],'linewidth',5)
    end
    %print(gcf,'-r300','-dpng','testPSD.png');
end



end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Low level function  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [SAfit SAforce fcorner] = spectral_analysis(fs,mean_z,x,plotflag)
% Calibrates the forces by analyzing the spectrum of a single trace in the
% manner described by te Velthuis et al. Bioph J 2010.
% This entails an iterative correction for aliasing and blurring, while
% fitting the integral of the PSD to yield the two free parameters (force and Rbead, if you wish)
%
% Input:
% - fs: measurement frequency in Hz
% - mean_z: mean extension in um
% - x: data trace in um
% - plotflag: to plot set 1 - to not plot set 0
%
%
% Output:
% - SAfit: (1) alpha in pN s/nm (2) kappa in pN/nm (3) good/bad (1/0)
% - SAforce: force according to ML PSD fit in pN
% - fcorner: corner frequency (according to atan fit)
%
%
% General remark: all analysis is performed in terms of frequency in Hz,
% the employed units are pN,nm and s/Hz everywhere



%%% PREPARE THE DATA %%%

%%% Get an estimate of the force and the corner frequency
kT = 4; %pN um
F_est = kT/1000*mean_z/std(x-mean(x))^2; %pN
f_c = calc_fcorner(F_est,mean_z); %Hz
fitgood = f_c < fs/2;

%%% Add a line to x to make the length 2^integer
x(end+1) = x(1);
%%% Subtract mean, work in units of nm
x = 1000*(x - mean(x));
mean_z = 1000 * mean_z;
N = length(x);

%%% Multipy by a hanning window - factor sqt(8/3) is required to maintain
%%% the same total power (the Matlab function 'hann' apparently doesn't care)
Hannwindow = hann(N);
x = x.*Hannwindow*sqrt(8/3);





%%% COMPUTE THE POWER SPECTRUM (PSD) %%%

[f PSD c] = calc_powersp(x,fs);
Nf = length(f);

%%% Save the original PSD for future reference
PSDold = PSD;



%%% ITERATIVELY FIND THE FORCE %%%
%%% Fit the integral of the power spectrum to get the model parameters.
%%% Next, perform the spectral correction for aliasing and blurring. Redo
%%% the fitting with the corrected spectrum. Repeat until the fitting
%%% parameters do not change anymore between iterations.

%%% NOTE: 'Hard-coded' maximum number of iterations
maxit = 20;

for loopind = 1:maxit
    
    %%% INTEGRATE AND FIT THE SPECTRUM %%%
    
    %%% Use trapezoidal numerical integration, note that first point is 0
    PSDintegral = cumtrapz(f,PSD);
    
    %%% Fit an arctangent to the integral of the PSD
    initguess = [max(PSDintegral) 1];
    [par res ef] = fminsearch(@(par) sum((PSDintegral-par(1)*atan(f/par(2))).^2),initguess);
    prefac = par(1); fcorner = par(2);
    arctanmodel = prefac*atan(f/fcorner);
    
    
    %%% From the prefac and the fcorner, calc alpha, kappa, force, Rbead
    alpha = kT/(prefac*2*pi^2*fcorner); %in pN s/nm
    kappa = 2*pi*fcorner*alpha; %in pN/nm
    force = kappa * mean_z; %in pN
    Rbead = alpha/(6*pi*10^-6); %Apparent bead size in mum assuming that
    % the viscosity of water, in Pa*s = N*s/m^2 is 0.001
    
    
    %%% Save the important parameters to see if the results converge
    forces(loopind) = force;
    alphas(loopind) = alpha;
    Rbeads(loopind) = Rbead;
    
    
    %%% If the difference with the previous iteration is negligible, stop
    %%% the for loop, we have our result. Otherwise, continue with the
    %%% spectral corrections. This difference is 'hard-coded'
    if loopind > 1
        if abs(force-forces(loopind-1)) < 1E-2 && abs(Rbead-Rbeads(loopind-1)) < 1E-2
            disp(['Converged after ' num2str(loopind) ' iterations'])
            break
        else
            if loopind == maxit
                disp(['SA method did not converge after ' num2str(maxit) ' iterations'])
                %%% NOTE: IF SO, PERHAPS AVERAGE OVER ALL PREVIOUS RESULTS?
            end
        end
    end
    
    warning off
    if plotflag == 1
        %%% Plot the atan fit of the integral of the PSD
        figure(3); clf; hold on; box on;
        plot(f,PSDintegral,'r','linewidth',2)
        plot(f,arctanmodel,'k--','linewidth',2)
        set(gca, 'fontsize', 14, 'linewidth', 1, 'fontweight', 'bold', 'TickLength',[0.02 0.02]);
        xlabel('Frequency (Hz)')
        ylabel('Integral of PSD (nm^2)')
        %set(gca,'yscale','log')
        %set(gca,'xscale','log')
        set(gca,'xlim',[f(1)/2 f(end)*1.1])
        title(['F_{var} = ' num2str(F_est,2) ' pN, F_{fit} = ' num2str(force,2) ' pN'])
        legend('Corrected PSD integral','Arctangent fit','Fancy model','location','east')
        legend(gca,'boxoff')
        %print(gcf,'-r300','-dpng','testSA1.png');
    end
    warning on
    
    %%% NEXT --- RUN ALL THE SPECTRAL CORRECTIONS %%%
    % see: te Velthuis et al. Bioph J 2010
    
    %%% Calculate the blurring window correction
    %%% NOTE: the sinc(x) function in matlab is sin(pi*x)/(pi*x)
    Cblur = sinc(f/fs).^2;
    
    %%% Calculate the aliasing corrections
    %%% The alias term is the product of the theoretical Fourier spectrum
    %%% and the blurring window, here for the spectra around fs and -fs
    ns = [-1 1];
    Calias = 0;
    for n = ns;
        fnow = f + n*fs;
        Lorentzian = kT/(2*pi^2*alpha)./(fnow.^2+fcorner^2);
        Pwindow = sinc(fnow/fs).^2;
        Calias = Calias + Lorentzian.*Pwindow;
    end
    
    %%% Correct the measured power spectrum by first subtracting the aliased
    %%% terms around fs and -fs and next by dividing through Cblur to account
    %%% for the blurring effect
    PSDcorr = (PSDold - Calias)./Cblur;
    PSD = PSDcorr;
    
    %%% NOTE: THIS INEVITABLY LEADS TO NEGATIVE POWER AT SOME FREQUENCIES!
    %%% Now set off the warning option in Matlab to keep the command window
    %%% clean
    warning off
    
    if plotflag == 1
        %%% Plot the uncorr and corrected PSD with Lorentzian model
        figure(4); clf; hold on; box on;
        plot(f,PSDold,'r','linewidth',1)
        plot(f,PSDcorr,'b','linewidth',1)
        plot(f,kT/(2*pi^2*alpha)./(f.^2+fcorner^2),'r-','linewidth',2)
        legend('Uncorrected PSD','Corrected PSD','Lorentzian model from integral fit','location','southwest')
        legend(gca,'boxoff')
        set(gca, 'fontsize', 14, 'linewidth', 1, 'fontweight', 'bold', 'TickLength',[0.02 0.02]);
        xlabel('Frequency (Hz)')
        ylabel('PSD (nm^2/Hz)')
        set(gca,'yscale','log')
        set(gca,'xscale','log')
        set(gca,'xlim',[f(1)/1.2 f(end)*1.2])
        title(['F_{var} = ' num2str(F_est,2) ' pN, F_{fit} = ' num2str(force,2) ' pN'])
        %print(gcf,'-r300','-dpng','testSA2.png');
    end
    
    warning on
end

%%% Display the iteration results if wanted
F_Rbead = [forces' Rbeads'];

%%% Save the results
SAfit = [alpha kappa fitgood];
SAforce = force;


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Low level function  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f,P,T] = calc_powersp(X,sampling_f)
% function [f,P,T] = calc_powersp(X,sampling_f)
% Calculates the powerspectrum P as function of frequency f
% of imput data X sampled at frequency sampling_f

global T

fNyq    =   sampling_f / 2;
delta_t =   1 / sampling_f;

time    =   [0 : delta_t : (length(X)-1)*delta_t]';
T       =   max(time);
f       =   ([1 : length(X)] / T)';

FT      =   delta_t*fft(X);
P       =   FT .* conj(FT) / T;

ind     =   find(f <= fNyq); % only to the Nyquist f
f       =   f(ind);
P       =   P(ind);


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Low level function  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function f_c = calc_fcorner(forces,mean_z)
% Runs a simple calculation to estimate the corner frequency in MT
% Hard-coded are the viscosity and bead-radius
% Input are the forces (vector in pN) and mean_z (vector in um)
% Output is the theoretical corner frequency f_c (vector in Hz)

eta      = 0.001;      %%% Viscosity of water, in Pa * s = N * s / m^2
%R_bead   = 0.5*10^(-6);    %%% Bead radius in m - MYONE
R_bead   = 1.4*10^(-6);    %%% Bead radius in m - M270

Length = mean_z*10^-6; %%% Tether extension in m
forces = forces*10^-12; %%% forces in N



f_c = 1/(2*pi).* forces ./ Length ./(6*pi*eta*R_bead);


%%% If necessary, plot the results
if 0
    plot(forces*10^12, f_c, 'bo')
    line([min(forces*10^12) max(forces*10^12)], [30 30], 'color', [1 0 0])
    set(gca, 'fontsize', 22, 'linewidth', 1, 'fontweight', 'bold', 'TickLength',[0.02 0.02]);
    xlabel('Force (pN)')
    ylabel('Corner frequency (Hz)')
    set(gca, 'yscale', 'log')
    set(gca, 'xscale', 'log')
end



end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Low level function  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function costfunc = costfunction_PSDfit(f,alpha,kappa,fs,PSD)
% Costfunction for PSD Max Likelihood fit according to
% the model by Lansdorp and Saleh (RSI 2012)

eta = 1; % As far as I know, this is constant anyway, so why calc?
PSDmodel = analytical_PSD_overdamped_bead(alpha,kappa,fs,f);

costfunc = eta .* (PSD./PSDmodel + log(PSDmodel));
costfunc = sum(costfunc);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Low level function  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function costfunc = costfunction_AVfit(tend,taus,alpha,kappa,AV)
% Costfunction for Allan variance Max Likelihood fit according to
% the model by Lansdorp and Saleh (RSI 2012)

eta = 1/2*(tend./taus-1);
AVmodel = analytical_AV_overdamped_bead(alpha,kappa,taus);

costfunc = eta .* (AV./AVmodel + log(AVmodel));
costfunc = sum(costfunc);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Low level function  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [Lp, Lc, Ffit, ext_fit] = inExtWLCfit(ext, forceIn, UseErr)
%
% Take a;ready determined forceIns and tether extensions
% and fit to the WLC equation to determine
% the persistence (Lp) and contour length (Lc)
%
% Arguments:
% fmax (optional) - forceIn values above this value will NOT be used in the WLC fit
% fmin (optional) - forceIn values below this value will NOT be used in the WLC fit
% UseErr (optional) - Set to "0" to not use errors; for a positive number,
% assume that level of relative error on the forceIn
%

if nargin < 2
    error('ERROR: Need to specify inputs!');
end
if nargin < 4
    disp('No error level specified - use no error bars');
    UseErr = 0;
end


%%% for the fits, select points based on fmax
forceIn_raw = forceIn;
ext_raw = ext;


control.fmax = 5;
control.fmin = 0;

ind = find(forceIn < control.fmax & forceIn > control.fmin);
forceIn = forceIn(ind);
ext = ext(ind);

if UseErr > 0
    
    %%%--- Fit the WLC model WITH ERRORS ---%%%
    %%% Use the model by Bouchiat, et al. Biophys J 76:409 (1999) %%%
    Lp0 = 0.045;    % Initial value for persistence length
    Lc0 = max(ext);   % Initial value for contour length
    
    forceIn_err = forceIn*UseErr;
    [par res ef] = fminsearch(@(par) sum( ((forceIn - wlc_7param(ext, par(1), par(2)))./forceIn_err).^2 )  , [Lp0 Lc0]);
    
    Lp = par(1);
    Lc = par(2);
    ext_fit = min(ext):0.01:max(ext);
    
    if ef == 1
        Ffit = wlc_7param(ext_fit, Lp, Lc);
    else
        Ffit = zeros(size(ext_fit));
        disp('ERROR: Fit did not converge!');
    end
    
    
else
    
    %%%--- Fit the WLC model WITHOUT ERRORS ---%%%
    %%% Use the model by Bouchiat, et al. Biophys J 76:409 (1999) %%%
    Lp0 = 0.045;    % Initial value for persistence length
    Lc0 = max(ext);   % Initial value for contour length
    
    
    [par res ef] = fminsearch(@(par) norm(forceIn - wlc_7param(ext, par(1), par(2)),2), [Lp0 Lc0]);
    Lp = par(1);
    Lc = par(2);
    ext_fit = min(ext):0.01:max(ext);
    
    if ef == 1
        Ffit = wlc_7param(ext_fit, Lp, Lc);
    else
        Ffit = zeros(size(ext_fit));
        disp('ERROR: Fit did not converge!');
    end
    
end


forceIn = forceIn_raw;
ext = ext_raw;


end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Low level function  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function s=exp2fit(t,f,caseval,lsq_val,options)

%  exp2fit solves the non-linear least squares problem exact
%  and using it as a start guess in a least square method
%  in cases with noise, of the specific exponential functions:
%  --- caseval = 1 ----
%  f=s1+s2*exp(-t/s3)
%
%  --- caseval = 2 (general case, two exponentials) ----
%  f=s1+s2*exp(-t/s3)+s4*exp(-t/s5)
%
%  --- caseval = 3 ----
%  f=s1*(1-exp(-t/s2)) %i.e., constraints between s1 and s2
%
%  Syntax: s=exp2fit(t,f,caseval) gives the parameters in the fitting
%  function specified by the choice of caseval (1,2,3).
%  t and f are (normally) vectors of the same size, containing
%  the data to be fitted.
%  s=exp2fit(t,f,caseval,lsq_val,options), using lsq_val='no' gives
%  the analytic solution, without least square approach (faster), where
%  options (optional or []) are produced by optimset, as used in lsqcurvefit.
%
%  This algorithm is using analytic formulas using multiple integrals.
%  Integral estimations are used as start guess in lsqcurvefit.
%  Note: For infinite lengths of t, and f, without noise
%  the result is exact.
%
% %--- Example 1:
% t=linspace(1,4,100)*1e-9;
% noise=0.02;
% f=0.1+2*exp(-t/3e-9)+noise*randn(size(t));
%
% %--- solve without startguess
% s=exp2fit(t,f,1)
%
% %--- plot and compare
% fun = @(s,t) s(1)+s(2)*exp(-t/s(3));
% tt=linspace(0,4*s(3),200);
% ff=fun(s,tt);
% figure(1), clf;plot(t,f,'.',tt,ff);
%
% %--- Example 2, Damped Harmonic oscillator:
% %--- Note: sin(x)=(exp(ix)-exp(-ix))/2i
% t=linspace(1,12,100)*1e-9;
% w=1e9;
% f=1+3*exp(-t/5e-9).*sin(w*(t-2e-9));
%
% %--- solve without startguess
% s=exp2fit(t,f,2,'no')
%
% %--- plot and compare
% fun = @(s,t) s(1)+s(2)*exp(-t/s(3))+s(4)*exp(-t/s(5));
% tt=linspace(0,20,200)*1e-9;
% ff=fun(s,tt);
% figure(1), clf;plot(t,f,'.',tt,real(ff));
% %--- evaluate parameters:
% sprintf(['f=1+3*exp(-t/5e-9).*sin(w*(t-2e-9))\n',...
% 'Frequency: w_fitted=',num2str(-imag(1/s(3)),3),' w_data=',num2str(w,3),'\n',...
% 'Damping: tau=',num2str(1/real(1/s(3)),3),'\n',...
% 'Offset: s1=',num2str(real(s(1)),3)])
%
%%% By Per Sundqvist january 2009.

[t,ix]=sort(t(:));%convert to column vector and sort
f=f(:);f=f(ix);

if nargin<4
    lsq_val='yes';%default, use lsq-fitting
end
if nargin<5
    options=optimset('TolX',1e-6,'TolFun',1e-8);%default
end
if nargin>=5
    if isempty(options)
        options=optimset('TolX',1e-6,'TolFun',1e-8);
    end
end
if length(t)<3
    error(['WARNING!', ...
        'To few data to give correct estimation of parameters!']);
end

%calculate help-variables
T=max(t)-min(t);t2=max(t);
tt=linspace(min(t),max(t),200);
ff=pchip(t,f,tt);
n=1;I1=trapz(tt,ff.*(t2-tt).^(n-1))/factorial(n-1);
n=2;I2=trapz(tt,ff.*(t2-tt).^(n-1))/factorial(n-1);
n=3;I3=trapz(tt,ff.*(t2-tt).^(n-1))/factorial(n-1);
n=4;I4=trapz(tt,ff.*(t2-tt).^(n-1))/factorial(n-1);

if caseval==1
    %--- estimate tau, s1,s2
    %--- Case: f=s1+s2*exp(-t/tau)
    tau=(12*I4-6*I3*T+I2*T^2)/(-12*I3+6*I2*T-I1*T^2);
    Q1=exp(-min(t)/tau);
    Q=exp(-T/tau);
    s1=2.*T.^(-1).*((1+Q).*T+2.*((-1)+Q).*tau).^(-1).*(I2.*((-1)+Q)+I1.* ...
        (T+((-1)+Q).*tau));
    s2=(2.*I2+(-1).*I1.*T).*tau.^(-1).*((1+Q).*T+2.*((-1)+Q).*tau).^(-1);
    s2=s2/Q1;
    sf0=[s1 s2 tau];
    fun = @(s,t) (s(1)*sf0(1))+(s(2)*sf0(2))*exp(-t/(s(3)*sf0(3)));
    s0=[1 1 1];
elseif caseval==3
    %--- estimate tau, s1
    %--- Case: f=s1*(1-exp(-t/tau))
    tau=(12*I4-6*I3*T+I2*T^2)/(-12*I3+6*I2*T-I1*T^2);
    s1=6.*T.^(-3).*((-2).*I3+I2.*(T+(-2).*tau)+I1.*T.*tau);
    sf0=[s1 tau];
    fun = @(s,t) (s(1)*sf0(1))*(1-exp(-t/(s(2)*sf0(2))));
    s0=[1 1];
elseif caseval==2
    %
    T=max(t)-min(t);t2=max(t);
    tt=linspace(min(t),max(t),200);
    ff=pchip(t,f,tt);
    n=1;J(n)=trapz(tt,ff.*(t2-tt).^(n-1))/factorial(n-1)/T^n;
    n=2;J(n)=trapz(tt,ff.*(t2-tt).^(n-1))/factorial(n-1)/T^n;
    n=3;J(n)=trapz(tt,ff.*(t2-tt).^(n-1))/factorial(n-1)/T^n;
    n=4;J(n)=trapz(tt,ff.*(t2-tt).^(n-1))/factorial(n-1)/T^n;
    n=5;J(n)=trapz(tt,ff.*(t2-tt).^(n-1))/factorial(n-1)/T^n;
    n=6;J(n)=trapz(tt,ff.*(t2-tt).^(n-1))/factorial(n-1)/T^n;
    n=7;J(n)=trapz(tt,ff.*(t2-tt).^(n-1))/factorial(n-1)/T^n;
    
    %
    p(1)=(1/2).*(J(2).^2+(-1).*J(1).*(J(3)+(-15).*(J(4)+(-6).*J(5)+14.*J(6) ...
        ))+(-15).*J(2).*(J(3)+2.*(J(4)+(-25).*J(5)+84.*J(6)))+120.*(J(3) ...
        .^2+J(3).*((-8).*J(4)+(-9).*J(5)+105.*J(6))+15.*(2.*J(4).^2+14.*J( ...
        5).^2+(-7).*J(4).*(J(5)+2.*J(6))))).^(-1).*(J(1).*(J(4)+(-15).*(J( ...
        5)+(-6).*J(6)+14.*J(7)))+(-1).*J(2).*(J(3)+(-120).*(J(5)+(-8).*J( ...
        6)+21.*J(7)))+15.*(J(3).^2+(-2).*J(3).*(7.*J(4)+(-7).*J(5)+(-120) ...
        .*J(6)+420.*J(7))+8.*(8.*J(4).^2+105.*J(5).*(J(5)+(-2).*J(6))+J(4) ...
        .*((-51).*J(5)+210.*J(7))))+sqrt(4.*(J(2).^2+(-1).*J(1).*(J(3)+( ...
        -15).*(J(4)+(-6).*J(5)+14.*J(6)))+(-15).*J(2).*(J(3)+2.*(J(4)+( ...
        -25).*J(5)+84.*J(6)))+120.*(J(3).^2+J(3).*((-8).*J(4)+(-9).*J(5)+ ...
        105.*J(6))+15.*(2.*J(4).^2+14.*J(5).^2+(-7).*J(4).*(J(5)+2.*J(6))) ...
        )).*((-1).*J(3).^2+J(2).*(J(4)+(-15).*(J(5)+(-6).*J(6)+14.*J(7)))+ ...
        15.*J(3).*(J(4)+2.*(J(5)+(-25).*J(6)+84.*J(7)))+(-120).*(J(4).^2+ ...
        J(4).*((-8).*J(5)+(-9).*J(6)+105.*J(7))+15.*(2.*J(5).^2+14.*J(6) ...
        .^2+(-7).*J(5).*(J(6)+2.*J(7)))))+(J(1).*(J(4)+(-15).*(J(5)+(-6).* ...
        J(6)+14.*J(7)))+(-1).*J(2).*(J(3)+(-120).*(J(5)+(-8).*J(6)+21.*J( ...
        7)))+15.*(J(3).^2+(-2).*J(3).*(7.*J(4)+(-7).*J(5)+(-120).*J(6)+ ...
        420.*J(7))+8.*(8.*J(4).^2+105.*J(5).*(J(5)+(-2).*J(6))+J(4).*(( ...
        -51).*J(5)+210.*J(7))))).^2));
    %
    p(2)=(-1/2).*(J(2).^2+(-1).*J(1).*(J(3)+(-15).*(J(4)+(-6).*J(5)+14.*J( ...
        6)))+(-15).*J(2).*(J(3)+2.*(J(4)+(-25).*J(5)+84.*J(6)))+120.*(J(3) ...
        .^2+J(3).*((-8).*J(4)+(-9).*J(5)+105.*J(6))+15.*(2.*J(4).^2+14.*J( ...
        5).^2+(-7).*J(4).*(J(5)+2.*J(6))))).^(-1).*((-1).*J(1).*(J(4)+( ...
        -15).*(J(5)+(-6).*J(6)+14.*J(7)))+J(2).*(J(3)+(-120).*(J(5)+(-8).* ...
        J(6)+21.*J(7)))+(-15).*(J(3).^2+(-2).*J(3).*(7.*J(4)+(-7).*J(5)+( ...
        -120).*J(6)+420.*J(7))+8.*(8.*J(4).^2+105.*J(5).*(J(5)+(-2).*J(6)) ...
        +J(4).*((-51).*J(5)+210.*J(7))))+sqrt(4.*(J(2).^2+(-1).*J(1).*(J( ...
        3)+(-15).*(J(4)+(-6).*J(5)+14.*J(6)))+(-15).*J(2).*(J(3)+2.*(J(4)+ ...
        (-25).*J(5)+84.*J(6)))+120.*(J(3).^2+J(3).*((-8).*J(4)+(-9).*J(5)+ ...
        105.*J(6))+15.*(2.*J(4).^2+14.*J(5).^2+(-7).*J(4).*(J(5)+2.*J(6))) ...
        )).*((-1).*J(3).^2+J(2).*(J(4)+(-15).*(J(5)+(-6).*J(6)+14.*J(7)))+ ...
        15.*J(3).*(J(4)+2.*(J(5)+(-25).*J(6)+84.*J(7)))+(-120).*(J(4).^2+ ...
        J(4).*((-8).*J(5)+(-9).*J(6)+105.*J(7))+15.*(2.*J(5).^2+14.*J(6) ...
        .^2+(-7).*J(5).*(J(6)+2.*J(7)))))+(J(1).*(J(4)+(-15).*(J(5)+(-6).* ...
        J(6)+14.*J(7)))+(-1).*J(2).*(J(3)+(-120).*(J(5)+(-8).*J(6)+21.*J( ...
        7)))+15.*(J(3).^2+(-2).*J(3).*(7.*J(4)+(-7).*J(5)+(-120).*J(6)+ ...
        420.*J(7))+8.*(8.*J(4).^2+105.*J(5).*(J(5)+(-2).*J(6))+J(4).*(( ...
        -51).*J(5)+210.*J(7))))).^2));
    %
    s2=3.*p(1).^(-1).*(p(1)+(-1).*p(2)).^(-1).*((-1).*J(2).*p(1)+(-4).*( ...
        J(2).*p(1).^2.*(2+5.*p(1))+5.*J(5).*(1+6.*p(1).*(1+2.*p(1))))+(( ...
        -1).*J(1).*p(1).*(1+4.*p(1).*(2+5.*p(1)))+J(2).*((-1)+12.*p(1) ...
        .^2.*(3+10.*p(1)))).*p(2)+J(3).*((-1)+8.*p(2)+12.*p(1).*(p(1).*(3+ ...
        p(1).*(10+(-20).*p(2)))+3.*p(2)))+(-4).*J(4).*((-2)+5.*p(2)+3.*p( ...
        1).*((-3)+10.*p(2)+20.*p(1).*(p(1)+p(2)))));
    %
    s3=3.*(p(1)+(-1).*p(2)).^(-1).*p(2).^(-1).*(J(3)+(-8).*J(4)+20.*J(5)+ ...
        J(2).*p(1)+(-8).*J(3).*p(1)+20.*J(4).*p(1)+(J(2)+(-36).*J(4)+120.* ...
        J(5)+(J(1)+(-36).*J(3)+120.*J(4)).*p(1)).*p(2)+(-4).*(9.*J(3)+( ...
        -60).*J(5)+(-2).*(J(1)+30.*J(4)).*p(1)+J(2).*((-2)+9.*p(1))).*p(2) ...
        .^2+20.*(J(2)+(-6).*J(3)+12.*J(4)+J(1).*p(1)+(-6).*(J(2)+(-2).*J( ...
        3)).*p(1)).*p(2).^3);
    %
    s1=6.*((-1)+(-3).*p(2)+(-1).*p(1).*(3+6.*p(2))).^(-1).*((-1).*J(3)+( ...
        1/2).*((s3+2.*s3.*p(1)+(-2).*J(1).*p(1)).*p(2)+(-2).*J(2).*(p(1)+ ...
        p(2))+s2.*p(1).*(1+2.*p(2))));
    %
    tau1=p(1)*T;
    tau2=p(2)*T;
    Q1=exp(-min(t)/tau1);
    Q2=exp(-min(t)/tau2);
    s2=s2/Q1;
    s3=s3/Q2;
    %
    sf0=[s1 s2 tau1 s3 tau2];
    fun = @(s,t) (s(1)*sf0(1))+...
        (s(2)*sf0(2))*exp(-t/(s(3)*sf0(3)))+...
        (s(4)*sf0(4))*exp(-t/(s(5)*sf0(5)));
    s0=[1 1 1 1 1];
end

%--- use lsqcurvefit if not lsq_val='no'
if isequal(lsq_val,'no')
    s=sf0;
else
    cond=1;
    while cond
        [s,RESNORM,RESIDUAL,EXIT]=lsqcurvefit(fun,s0,t,f,[],[],options);
        cond=not(not(EXIT==0));
        s0=s;
    end
    s=s0.*sf0;
end


end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Low level function  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Fwlc] = wlc_7param(ext, Lp, Lc)
%
% Given a vector of extensions and the parameter
% Lp = persistence length  (in micro-m)
% Lc = contour length      (in micro-m)
% T  = absolute temperature (in Kelvin)
% This function returns the forces computed from a
% 7 parameter model of the WLC,
% using the model by Bouchiat, et al. Biophys J 76:409 (1999)

T = 23+273.15;
kT = 1.3806503 *10^(-23)*T*10^18; % k_B T in units of pN micro-m
z_scaled = ext/Lc;

a  = zeros(6,1);
a(1) = 1;
a(2) = -0.5164228;
a(3) = -2.737418;
a(4) = 16.07497;
a(5) = -38.87607;
a(6) = 39.49944;
a(7) = -14.17718;


Fwlc = 1./(4*(1 - z_scaled).^2)  - 1/4;
for i=1:7
    Fwlc = Fwlc + a(i) .* z_scaled.^(i);
end

Fwlc = Fwlc * kT/Lp;

end