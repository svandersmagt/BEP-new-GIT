function bead = subtractZOffsetHighestPlateau(zmags,plat,bead,configVariable) 
%%% This function is called when no seperate measurement for determining
%%% the z-offset was supplied, and the z-offset was not yet subtracted. It
%%% uses the plateau with the weakest magnet force to determine the
%%% z-offset, in the same way step 1, zOffsetCalculation.m does.

%%% Input: (zmags,plat,bead,configVariable)
%%% - vector containing the magnetheight for every plateau
%%% - struct containing the first and last index for each plateau
%%% - struct containing the time trace and position traces of the beads

%%% Output: bead
%%% - returns bead, with z-offset subtracted
%%
    plotThings = configVariable.plotThings;
    [~,highest] = max(zmags);
    highestPlateauZ = bead(1).z(plat(highest).first:plat(highest).last);
    highestPlateauTime = bead(1).time(plat(highest).first:plat(highest).last);
    zOffsets = [];
    nSmooth=100;

    for i=1 %there is only one bead in this example

        %%% Smooth and find minimum
        smoothZ = smooth(highestPlateauZ, nSmooth, 'moving');
        [minZ, ind] = min(smoothZ);
        zOffsets = [zOffsets minZ];
        bead(i).z = bead(i).z - minZ;
        display('Z-offset subtracted using highest plateau');
        
        if plotThings;
            figure(1); clf; hold on; box on;
            plot(highestPlateauTime, highestPlateauZ, 'k-', 'linewidth', 1);
            plot(highestPlateauTime, smoothZ, 'r-', 'linewidth', 2);
            plot(highestPlateauTime(ind), smoothZ(ind), 'bx', 'linewidth', 2, 'markersize', 20);
            xlabel('Time (s)'); ylabel('z (nm)');
            title(['Z-offset finding for bead # ' num2str(i)]);
            legend('z','Smoothed z','Lowest point');
        end
    end
end