function [bead] = analyzeData(bead, nPlat, plat, configVariable)
%%% STEP 4: ANALYZE THE BEAD DATA FOR EACH PLATEAU

%%% Some variables ----------------------------------------------
kT = configVariable.kT; %pN nm
viscosity = configVariable.viscosity; %viscosity in pN s/nm^2
beadRadius = configVariable.beadRadius; %bead radius in nm
nBlock = configVariable.nBlock; %number of blocks to use for blocked fourier transform
sampleFreq = configVariable.sampleFreq;
plotThings = configVariable.plotThings;

for i=1;

    bead(i).ext = zeros(1,nPlat);
    bead(i).realTimeForceX = zeros(1,nPlat);
    bead(i).realTimeForceY = zeros(1,nPlat);
    bead(i).ZhongboYForce = zeros(1,nPlat);
    bead(i).DaldropXForce = zeros(1,nPlat);
    bead(i).DaldropXRadius = zeros(1,nPlat);
    bead(i).DaldropYForce = zeros(1,nPlat);
    bead(i).DaldropYRadius = zeros(1,nPlat);
    bead(i).ZhongboDaldropForceDifference = zeros(1,nPlat);
    bead(i).DaldropXYForceDifference = zeros(1,nPlat);
    bead(i).DaldropXYRadiusDifference = zeros(1,nPlat);


    %%% Use the script "analyzeTrace" to determine the force for each trace
    for k=1:nPlat

        [ext, realTimeForceX, realTimeForceY, ~, ZhongboYForce, ~, DaldropXForce, DaldropXRadius, ~, DaldropYForce, DaldropYRadius]=...
        analyzeTrace(...
            bead(i).x(plat(k).first:plat(k).last),... 
            bead(i).y(plat(k).first:plat(k).last),...
            bead(i).z(plat(k).first:plat(k).last),...
        sampleFreq, beadRadius, kT, viscosity, nBlock, plotThings);

        bead(i).ext(k) = ext;
        bead(i).realTimeForceX(k) = realTimeForceX;
        bead(i).realTimeForceY(k) = realTimeForceY;
        bead(i).ZhongboYForce(k) = ZhongboYForce;
        bead(i).DaldropXForce(k) = DaldropXForce;
        bead(i).DaldropXRadius(k) = DaldropXRadius;
        bead(i).DaldropYForce(k) = DaldropYForce;
        bead(i).DaldropYRadius(k) = DaldropYRadius;
        bead(i).ZhongboDaldropForceDifference(k) = abs(DaldropYForce - ZhongboYForce);
        bead(i).DaldropXYForceDifference(k) = abs(DaldropYForce - DaldropXForce);
        bead(i).DaldropXYRadiusDifference(k) = abs(DaldropYRadius - DaldropXRadius);

        display(['Finished working on bead # ' num2str(i) ', plateau number ' num2str(k)])
    end
end