% % % This script was used to be able to tweak the algorithm and test the 
% % % outcome against the values produced by Daldrops Labview implementation.
% % % It uses the data provided by Daldrop.
%%
clear all; clc; close all;

%%% Datafile and aquisition frequency
tracesFile = '..\BEPmaarnietgit\DaldropData.txt';
sampleFreq = 2800; %%% Acquisition frequency in Hz

%%% Constants (pN and nm) 
kT = 4.1; %pN nm
viscosity = 10E-10; %Viscosity in pN s/nm^2
beadRadius = 515; %Bead radius in nm
nBlock = 1; %Number of blocks used for blocked powerspectrum
plot = true;

%%% Read in data
data = load(tracesFile);    

%%% Parse the bead data
for i=1; 
    bead(i).time = 1:length(data(:,1));
    bead(i).x = data(:,1)*1000;
    bead(i).y = data(:,2)*1000;
    bead(i).z = data(:,3)*1000; %data is in mum, but this script uses nm
end

for i=1;
    %%% Analyse data in 3 ways: 
    %%% -Zhongbo's method of fitting the spectrum to the analytical solution in the y-direction.
    %%% -Daldrop's method of fitting the spectrum to the analytical solution in the y-direction.
    %%% -Daldrop's method of fitting the spectrum to the analytical solution in the x-direction, taking bead rotations into account.
    %%% These three methods are applied in analyzeTrace()
    [ext, realTimeForceX, realTimeForceY, ZhongboYFit, ZhongboYForce,...
        DaldropXFit, DaldropXForce, DaldropXRadius, DaldropYFit,...
        DaldropYForce, DaldropYRadius]= analyzeTrace(bead(i).x,... 
            bead(i).y,...
            bead(i).z,...
            sampleFreq, beadRadius, kT, viscosity, nBlock, plot);

    bead(i).ext = ext;
    bead(i).realTimeForceX = realTimeForceX;
    bead(i).realTimeForceY = realTimeForceY;
    bead(i).ZhongboYForce = ZhongboYForce;
    bead(i).DaldropXForce = DaldropXForce;
    bead(i).DaldropXRadius = DaldropXRadius;
    bead(i).DaldropYForce = DaldropYForce;
    bead(i).DaldropYRadius = DaldropYRadius;

    %%% Display results for easy comparison
    display('least squares fit:')
    display(['Zhongbo method: F = ' num2str(bead(1).ZhongboYForce) 'pN, R = ' num2str(beadRadius) 'nm'])
    display(['Daldrop method y-direction: F = ' num2str(bead(1).DaldropYForce) 'pN, R = ' num2str(bead(1).DaldropYRadius) 'nm'])
    display(['Daldrop method x-direction: F = ' num2str(bead(1).DaldropXForce) 'pN, R = ' num2str(bead(1).DaldropXRadius) 'nm'])
    display(['Labview implementation of Daldrop y-direction: F = 4.9115pN, R = 512nm'])
    display(['Labview implementation of Daldrop x-direction: F = 4.9223pN, R = 514nm'])
end