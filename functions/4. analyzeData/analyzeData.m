function [bead, forcesExponentialFit] = analyzeData(bead, nPlat, plat, zmags, configVariable)
%%% STEP 4: ANALYZE THE BEAD DATA FOR EACH PLATEAU
%%% --------------------------------------------------------------------
%%% In this step the bead traces are analyzed. See analyzeTrace.m and
%%% fitExponentialToForces.m

%%% Input: (bead, nPlat, plat, zmags, configVariable)
%%% - struct containing the time trace and position traces of the beads
%%% - number of plateaus
%%% - struct containing the first and last index for each plateau
%%% - vector containing the magnetheight for every plateau

%%% Output: [bead, forcesExponentialFit]
%%% - result of the fits is added to the bead struct and returned
%%% - struct containing the result of the double exponential fit to the
%%%   forces vs. magnetheight
    %%
    %%% Some variables from the configuration file
    kT = configVariable.kT; %pN nm
    viscosity = configVariable.viscosity; %viscosity in pN s/nm^2
    beadRadius = configVariable.beadRadius; %bead radius in nm
    nBlock = configVariable.nBlock; %number of blocks to use for blocked fourier transform
    sampleFreq = configVariable.sampleFreq;
    plotThings = configVariable.plotThings;

    for i=1;

        bead(i).ext = zeros(1,nPlat);
        bead(i).realTimeForceLong = zeros(1,nPlat);
        bead(i).realTimeForceShort = zeros(1,nPlat);
        bead(i).ZhongboForceShort = zeros(1,nPlat);
        bead(i).DaldropForceLong = zeros(1,nPlat);
        bead(i).DaldropRadiusLong = zeros(1,nPlat);
        bead(i).DaldropForceShort = zeros(1,nPlat);
        bead(i).DaldropRadiusShort = zeros(1,nPlat);
        bead(i).ZhongboDaldropForceDifference = zeros(1,nPlat);
        bead(i).DaldropXYForceDifference = zeros(1,nPlat);
        bead(i).DaldropXYRadiusDifference = zeros(1,nPlat);


        %%% Use the script "analyzeTrace" to determine the force for each trace
        for k=1:nPlat

            [ext, realTimeForceLong, realTimeForceShort, ~, ZhongboForceShort, ~, ...
                DaldropForceLong, DaldropRadiusLong, ~, DaldropForceShort, DaldropRadiusShort]=...
            analyzeTrace(...
                bead(i).long(plat(k).first:plat(k).last),... 
                bead(i).short(plat(k).first:plat(k).last),...
                bead(i).z(plat(k).first:plat(k).last),...
            sampleFreq, beadRadius, kT, viscosity, nBlock, plotThings);

            bead(i).ext(k) = ext;
            bead(i).realTimeForceLong(k) = realTimeForceLong;
            bead(i).realTimeForceShort(k) = realTimeForceShort;
            bead(i).ZhongboForceShort(k) = ZhongboForceShort;
            bead(i).DaldropForceLong(k) = DaldropForceLong;
            bead(i).DaldropRadiusLong(k) = DaldropRadiusLong;
            bead(i).DaldropForceShort(k) = DaldropForceShort;
            bead(i).DaldropRadiusShort(k) = DaldropRadiusShort;
            bead(i).ZhongboDaldropForceDifference(k) = abs(DaldropForceShort - ZhongboForceShort);
            bead(i).DaldropXYForceDifference(k) = abs(DaldropForceShort - DaldropForceLong);
            bead(i).DaldropXYRadiusDifference(k) = abs(DaldropRadiusShort - DaldropRadiusLong);

            display(['Finished working on bead # ' num2str(i) ', plateau number ' num2str(k)]);
        end
        display('Fit force vs. magnet height to double exponential');
        [forcesExponentialFit] = fitExponentialToForces(bead, zmags, configVariable);
    end
end