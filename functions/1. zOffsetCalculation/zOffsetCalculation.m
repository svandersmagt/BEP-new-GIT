function [] = zOffsetCalculation(configVariable)
%%% STEP 1: FIND OFFSETS IN Z-DIRECTION PER BEAD  
%%% ---------------------------------------------------------------
%%% This section will be skipped if it is indicated in the configuration
%%% file that z-offsets were already saved to a file or already subtracted
%%% from the z-trace.

%%% Input: (configVariable)

%%% Output: []
%%% Saves z-offset data per bead to a file.
%%
    plotThings = configVariable.plotThings;

    %%% Write here some descriptive text about your data
    %%% ---
    %%% Example data set on FX measurements; M270 beads, 21 kbp DNA
    tracesFile = configVariable.zOffsetDataFile;
    outputFile = configVariable.zOffsetOutputFile;

    %%% Read in and parse bead data
    data = load(tracesFile);

    for i=1 %there is only one bead in this example
        if configVariable.firstColumnIsTime;
            bead(i).time = 1:length(data(:,1));
            bead(i).z = data(:,4);
        else
            bead(i).time = 1:length(data(:,1));
            bead(i).z = data(:,3);
    end

    zOffsets = [];
    nSmooth=100;

    for i=1 %there is only one bead in this example

        %%% Smooth and find minimum
        smoothZ = smooth(bead(i).z, nSmooth, 'moving');
        [minZ, ind] = min(smoothZ);
        zOffsets = [zOffsets minZ];

        if plotThings;
            figure(1); clf; hold on; box on;
            plot(bead(i).time, bead(i).z, 'k-', 'linewidth', 1);
            plot(bead(i).time, smoothZ, 'r-', 'linewidth', 2);
            plot(bead(i).time(ind), smoothZ(ind), 'bx', 'linewidth', 2, 'markersize', 20);
            xlabel('Time (s)'); ylabel('z (um)');
            title(['Z-offset finding for bead # ' num2str(i)]);
            legend('z','Smoothed z','Lowest point');
        end
    end

    %%% Save the data
    display('Save offsets to file')
    foo = [zOffsets'];
    save(outputFile, 'foo', '-ascii')
end