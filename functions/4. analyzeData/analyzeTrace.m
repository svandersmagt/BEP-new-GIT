function [ext, realTimeForceX, realTimeForceY, ZhongboYFit, ZhongboYForce, DaldropXFit, DaldropXForce, DaldropXRadius, DaldropYFit, DaldropYForce, DaldropYRadius]...
    = analyzeTrace(x, y, z, fs, R, kT, eta, nBlock, plot)
%%% Function to analyze magnetic tweezers time traces in three different ways:
%%% -Zhongbo's method of fitting the spectrum to the analytical solution in the y-direction.
%%% -Daldrop's method of fitting the spectrum to the analytical solution in the y-direction.
%%% -Daldrop's method of fitting the spectrum to the analytical solution in the x-direction, taking bead rotations into account.

%%% Input: (x, y, z, fs, R, kT, eta, nBlock)
%%% - trace x in nm
%%% - trace y in nm
%%% - trace z in nm
%%% - fs in Hz, sampling frequency
%%% - R in nm, initial guess for the bead Radius
%%% - kT in pN nm
%%% - eta, viscosity in pN s/nm^2
%%% - nBlock, number of blocks used for blocked powerspectrum

%%% Output: [ext, FxReal, FyReal, ZhongboYFit, ZhongboYForce, DaldropXFit, DaldropXForce, DaldropXRadius, DaldropYFit, DaldropYForce, DaldropYRadius]
%%% - results op analysis in real time, using Zhongbo's method, and using Daldrop's methods
%%

%%% Analyze the real time fluctuations
meanX = mean(x);
meanY = mean(y);
meanZ = mean(z);
stdX  = std(x);
stdY  = std(y);

realTimeForceX = kT*meanZ./stdX^2;
realTimeForceY = kT*meanZ./stdY^2;
ext = mean(sqrt((x-meanX).^2 + (y-meanY).^2 + z.^2));

%%% The three methods are implemented in their own functions:
fprintf('Zhongbo Fit \n');
[ZhongboYFit, ZhongboYForce] = analyzeZhongboY(fs,ext,y-meanY,kT,R,eta,nBlock,plot);
fprintf('Daldrop X Fit \n');
[DaldropXFit, DaldropXForce, DaldropXRadius] = analyzeDaldropX(fs,ext,x-meanX,y-meanY,R,kT,eta,nBlock,plot);
fprintf('Daldrop Y Fit \n');
[DaldropYFit, DaldropYForce, DaldropYRadius] = analyzeDaldropY(fs,ext,y-meanY,R,kT,eta,nBlock,plot);

end