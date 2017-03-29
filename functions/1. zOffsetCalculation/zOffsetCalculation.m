function [] = zOffsetCalculation(configVariable)
%%% STEP 1: FIND OFFSETS IN Z-DIRECTION PER BEAD  
%%% ---------------------------------------------------------------
%%% skip this section if your data has offsets already subtracted and
%%% skip this section if the offsets belonging to your data are already
%%% stored in a file
%%% ---------------------------------------------------------------

plotThings = configVariable.plotThings;

%%% Write here some descriptive text about your data
%%% ---
%%% Example data set on FX measurements; M270 beads, 21 kbp DNA
tracesFile = configVariable.zOffsetDataFile;
outputFile = configVariable.zOffsetOutputFile;

%%% Read in and parse bead data
data = load(tracesFile);

for i=1 %there is only one bead in this example
    bead(i).time = 1:length(data(:,1));
    bead(i).y = data(:,1);
    bead(i).z = data(:,2);
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

%%% Plot the offsets
if plotThings;
    figure(2); clf; hold on; box on;
    plot(1, zOffsets, 'bo', 'linewidth', 2, 'markersize', 5);
    title('Z-offset per bead');
    xlabel('Bead #'); ylabel('z-offset (um)');
    legend('Offset');
end

%%% Save the data
display('Save offsets to file')
foo = [(1)' zOffsets'];
save(outputFile, 'foo', '-ascii')
end