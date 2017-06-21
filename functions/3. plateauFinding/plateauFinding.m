function [plat, zmags, nPlat] = plateauFinding(zmag, configVariable)
%%% STEP 3: FINDING PLATEAUS IN THE MOTOR DATA
%%% ---------------------------------------------------------------
%%% Figure out where the magnets are moving, from the motor file
%%% This determines the "plateaus", where the magnet height is
%%% constant and where we we want to analyze the forces

%%% Input: (zmag, configVariable)
%%% - vector containing motor data of the magnet height

%%% Output: [plat, zmags, nPlat]
%%% - struct containing the first and last index for each plateau
%%% - vector containing the magnetheight for every plateau
%%% - number of plateaus
%%
    plotThings = configVariable.plotThings;

    %%% Some values that seem to work fine
    nSmoothZmag = 200;
    nSmoothdZmag = 200;
    threshold = 10^(-5); %threshold to determine where it is moving
    nMinPointsPlat = 100; %minimum number of points in a plateau

    %%% Smooth the motor data
    zmagSmooth = smooth(zmag, nSmoothZmag, 'moving');
    diffZmag = smooth(diff(zmagSmooth), nSmoothdZmag, 'lowess');

    nFirst =1; nLast  =1; %initializing

    %%% Check whether we are starting in a plateau
    if (abs(diffZmag(1)) < threshold && abs(diffZmag(2)) < threshold)
        tPlat(1).first = 1;
        nFirst = nFirst +1;
    end

    %%% Find plateaus inbetween
    for i=2:length(diffZmag)
        if (abs(diffZmag(i)) < threshold && abs(diffZmag(i-1)) > threshold)
            tPlat(nFirst).first = i;
            nFirst = nFirst +1;
        end

        if (abs(diffZmag(i)) > threshold && abs(diffZmag(i-1)) < threshold)
            tPlat(nLast).last = i;
            nLast = nLast +1;
        end
    end

    %%% Check whether we are ending in a plateau
    if (abs(diffZmag(end)) < threshold && abs(diffZmag(end-1)) < threshold)
        tPlat(nLast).last = length(diffZmag);
    end

    %%% Plateaus that are too short are thrown away
    nPlat = nFirst-1; 
    count = 1;
    Ngoodplat = 0;
    for j = 1:nPlat
        if length(tPlat(j).first:tPlat(j).last) > nMinPointsPlat
            plat(count).first = tPlat(j).first;
            plat(count).last = tPlat(j).last;

            Ngoodplat = Ngoodplat + 1;
            count = count + 1;

        else
            display(['Plateau # ' num2str(j) ' has only ' ...
                num2str(length(nMinPointsPlat)) ' points!' ])
        end
    end

    nPlat = Ngoodplat; %number of useable plateaus
    display(['Found ' num2str(nPlat) ' Zmag plateaus']);

    %%% Save the magnet height information
    zmags = zeros(1,nPlat);
    for j = 1:nPlat
        zmags(j) = zmagSmooth(plat(j).last);
    end

    if plotThings;
        figure(2);clf; hold on; box on;
        platind = find(abs(diffZmag) < threshold);
        faketime = 1:length(zmagSmooth);
        plot(faketime , zmagSmooth, 'b-')
        plot(faketime(platind), zmagSmooth(platind), 'r.')
        for j = 1:nPlat
            plot(faketime(plat(j).first) , zmagSmooth(plat(j).first), 'ko', 'markersize', 5)
            plot(faketime(plat(j).last) , zmagSmooth(plat(j).last), 'mo', 'markersize', 5)
        end
        xlabel('Time (s)'); ylabel('Magnet height (mm)')
        title(['Result of plateau finding; Found ' num2str(nPlat) ' plateaus.' ])
        legend( 'Motor data','Plateaus','Plateau start','Plateau end')
    end
end
