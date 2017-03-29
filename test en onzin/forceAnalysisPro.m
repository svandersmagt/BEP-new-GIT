
    
%%
%%% STEP 1: FIND OFFSETS IN Z-DIRECTION PER BEAD  
%%% ---------------------------------------------------------------
%%% skip this section if your data has offsets already subtracted and
%%% skip this section if the offsets belonging to your data are already
%%% stored in a file
%%% ---------------------------------------------------------------

clear all; clc; close all; format compact;
plotThings = false;

%%% Write here some descriptive text about your data
%%% ---
%%% Example data set on FX measurements; M270 beads, 21 kbp DNA
traces_file = '..\BEPmaarnietgit\offset.txt';
output_name = '..\BEPmaarnietgit\FX_offsets.txt';

%%% Read in and parse bead data
data = load(traces_file);

for i=1 %there is only one bead in this example
    bead(i).time = 1:length(data(:,1));
    bead(i).y = data(:,1);
    bead(i).z = data(:,2);
end

z_offsets = [];
Nsmooth=100;

for i=1 %there is only one bead in this example

    %%% Smooth and find minimum
    smooth_z = smooth(bead(i).z, Nsmooth, 'moving');
    [min_z, ind] = min(smooth_z);
    z_offsets = [z_offsets min_z];

    if plotThings;
        figure(1); clf; hold on; box on;
        plot(bead(i).time, bead(i).z, 'k-', 'linewidth', 1);
        plot(bead(i).time, smooth_z, 'r-', 'linewidth', 2);
        plot(bead(i).time(ind), smooth_z(ind), 'bx', 'linewidth', 2, 'markersize', 20);
        set(gca, 'fontsize', 16, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02]);
        xlabel('Time (s)'); ylabel('z (um)');
        title(['Bead # ' num2str(i)]);
    end
end

%%% Plot the offsets
if plotThings;
    figure(2); clf; hold on; box on;
    plot(1, z_offsets, 'bo', 'linewidth', 2, 'markersize', 20)
    set(gca, 'fontsize', 16, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02]);
    xlabel('Bead #'); ylabel('z-offset (um)')
end

%%% Save the data
foo = [(1)' z_offsets'];
save(output_name, 'foo', '-ascii')

%%
%%% STEP 2: LOAD DATA AND SUBTRACT Z-OFFSET

clear all; clc; close all;
plotThings = false;
offsetZ = true; %set to false if your data has offset already subtracted

%%% Write here some descriptive text about your data
%%% ---
%%% Example of force extension data
%%% 21 kbp DNA, M270 beads, 1 mm gap verticaly oriented magnets
traces_file = '..\BEPmaarnietgit\bead.txt';
motors_file = '..\BEPmaarnietgit\bead_motors.txt';
zoffsets_file = '..\BEPmaarnietgit\FX_offsets.txt';
freq = 60; %acquisition frequency in Hz

%%% Read in data
data = load(traces_file);
zmag = load(motors_file);

for i=1;
    bead(i).time = 1:length(data(:,1));
    bead(i).x = data(:,1)*1000; %nm
    bead(i).y = data(:,1)*1000; %only one trace in this example, so x and y are the same
    bead(i).z = data(:,2)*1000;
end

if offsetZ;
    %%% Read in the previously determined z-offsets from file
    zoff_data = load(zoffsets_file);
    zoffsets = zoff_data(:,2)*1000; %nm
    
    %%% Subtract previously determined z-offsets
    if length(zoffsets) == length(bead);
        for i=1;
            bead(i).z = bead(i).z - zoffsets(i);
        end
        display('Subtracted pre-determined z-offsets.')
    else
        display(['Number of beads in the big data file (' num2str(length(bead)) ') does not agree with the number of beads in the zoffset file (' num2str(length(zoffsets)) ')'])
    end
end

%%% STEP 3: FINDING PLATEAUS IN THE MOTOR DATA
%%% ---------------------------------------------------------------
%%% Figure out where the magnets are moving, from the motor file
%%% This determines the "plateaus", where the magnet height is
%%% constant and where we we want to analyze the forces
%%% ---------------------------------------------------------------

clc;

%%% Some values that seem to work fine
Nsmooth_zmag = 200;
Nsmooth_dzmag = 200;
small = 10^(-5); %threshold to determine where it is moving
Nmin_points_plat = 100; %minimum number of points in a plateau

%%% Smooth the motor data
zmag_smooth = smooth(zmag, Nsmooth_zmag, 'moving');
diff_zmag = smooth(diff(zmag_smooth), Nsmooth_dzmag, 'lowess');

Nfirst =1; Nlast  =1; %initializing

%%% Check whether we are starting in a plateau
if (abs(diff_zmag(1)) < small && abs(diff_zmag(2)) < small)
    tplat(1).first = 1;
    Nfirst = Nfirst +1;
end

%%% Find plateaus inbetween
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

%%% Plateaus that are too short are thrown away
Nplat = Nfirst-1; 
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

Nplat = Ngoodplat; %number of useable plateaus
display(['Found ' num2str(Nplat) ' Zmag plateaus']);

%%% Save the magnet height information
zmags = zeros(1,Nplat);
for j = 1:Nplat
    zmags(j) = zmag_smooth(plat(j).last);
end

if plotThings;
    figure(1);clf; hold on; box on;
    platind = find(abs(diff_zmag) < small);
    faketime = 1:length(zmag_smooth);
    plot(faketime , zmag_smooth, 'b-')
    plot(faketime(platind), zmag_smooth(platind), 'r.')
    for j = 1:Nplat
        plot(faketime(plat(j).first) , zmag_smooth(plat(j).first), 'ko', 'markersize', 10)
        plot(faketime(plat(j).last) , zmag_smooth(plat(j).last), 'mo', 'markersize', 10)
        
    end
    set(gca, 'fontsize', 14, 'fontweight', 'bold', 'linewidth', 1,'TickLength',[0.02 0.02])
    xlabel('Time'); ylabel('Zmag')
    title(['Red = plateaus; Black / magneta circles = start / stop; Found ' num2str(Nplat) ' plateaus.' ])
end

pause;

%%% STEP 4: ANALYZE THE BEAD DATA FOR EACH PLATEAU

%%% Some variables ----------------------------------------------
kT = 4.1; %pN nm
eta = 10E-10; %viscosity in pN s/nm^2
Rbead = 1400; %bead radius in nm
nBlock = 5; %number of blocks to use for blocked fourier transform

for i=1;

    bead(i).ext = zeros(1,Nplat);
    bead(i).FxReal = zeros(1,Nplat);
    bead(i).FyReal = zeros(1,Nplat);
    bead(i).ZhongboYForce = zeros(1,Nplat);
    bead(i).DaldropXForce = zeros(1,Nplat);
    bead(i).DaldropXRadius = zeros(1,Nplat);
    bead(i).DaldropYForce = zeros(1,Nplat);
    bead(i).DaldropYRadius = zeros(1,Nplat);


    %%% Use the script "analyzeTrace" to determine the force for each trace
    for k=1:Nplat

        [ext, FxReal, FyReal, ZhongboYFit, ZhongboYForce, DaldropXFit, DaldropXForce, DaldropXRadius, DaldropYFit, DaldropYForce, DaldropYRadius]=...
        analyzeTrace(...
            bead(i).x(plat(k).first:plat(k).last),... 
            bead(i).y(plat(k).first:plat(k).last),...
            bead(i).z(plat(k).first:plat(k).last),...
        freq, Rbead, kT, eta, nBlock, plotThings);

        bead(i).ext(k) = ext;
        bead(i).FxReal(k) = FxReal;
        bead(i).FyReal(k) = FyReal;
        bead(i).ZhongboYForce(k) = ZhongboYForce;
        bead(i).DaldropXForce(k) = DaldropXForce;
        bead(i).DaldropXRadius(k) = DaldropXRadius;
        bead(i).DaldropYForce(k) = DaldropYForce;
        bead(i).DaldropYRadius(k) = DaldropYRadius;

        display(['Finished working on bead # ' num2str(i) ', plateau number ' num2str(k)])
        %pause;
    end
end