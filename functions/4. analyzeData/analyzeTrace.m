function [extensionDNA, realTimeForceLong, realTimeForceShort, ZhongboFitShort,...
    ZhongboForceShort, DaldropFitLong, DaldropForceLong, DaldropRadiusLong,...
    DaldropFitShort, DaldropForceShort, DaldropRadiusShort]...
    = analyzeTrace(long, short, z, sampleFreq, beadRadius, kT, viscosity, nBlock, plotThings)
%%% Function to analyze magnetic tweezers time traces in three different ways:
%%% -Zhongbo's method of fitting the spectrum to the analytical solution in the short pendulum direction.
%%% -Daldrop's method of fitting the spectrum to the analytical solution in the short pendulum direction.
%%% -Daldrop's method of fitting the spectrum to the analytical solution in the long pendulum direction, taking bead rotations into account.

%%% Input: (long, short, z, sampleFreq, beadRadius, kT, viscosity, nBlock, plotThings)
%%% - trace in long pendulum direction in nm
%%% - trace in short pendulum direction in nm
%%% - trace z in nm
%%% - sampling frequency in Hz
%%% - initial guess for the bead Radius in nm
%%% - kT in pN nm
%%% - viscosity in pN s/nm^2
%%% - number of blocks used for blocked powerspectrum
%%% - show plots

%%% Output: [ext, realTimeForceLong, realTimeForceShort, ZhongboFitShort,...
%%% ZhongboForceShort, DaldropFitLong, DaldropForceLong, DaldropRadiusLong,...
%%% DaldropFitShort, DaldropForceShort, DaldropRadiusShort]
%%% - results op analysis in real time, using Zhongbo's method, and using Daldrop's methods
%%
    %%% Analyze the real time fluctuations
    meanLong = mean(long);
    meanShort = mean(short);
    meanZ = mean(z);
    stdLong  = std(long);
    stdShort  = std(short);

    realTimeForceLong = kT*meanZ./stdLong^2;
    realTimeForceShort = kT*meanZ./stdShort^2;
    extensionDNA = mean(sqrt((long-meanLong).^2 + (short-meanShort).^2 + z.^2));

    %%% The three methods are implemented in their own functions:
    fprintf('Zhongbo Short Pendulum Fit \n');
    [ZhongboFitShort, ZhongboForceShort] = analyzeZhongboShort(sampleFreq,extensionDNA,...
        short-meanShort,kT,beadRadius,viscosity,nBlock,plotThings);
    
    fprintf('Daldrop Long Pendulum Fit \n');
    [DaldropFitLong, DaldropForceLong, DaldropRadiusLong] = analyzeDaldropLong(sampleFreq,...
        extensionDNA,long-meanLong,short-meanShort,beadRadius,kT,viscosity,nBlock,plotThings);
    
    fprintf('Daldrop Short Pendulum Fit \n');
    [DaldropFitShort, DaldropForceShort, DaldropRadiusShort] = analyzeDaldropShort(...
        sampleFreq,extensionDNA,short-meanShort,beadRadius,kT,viscosity,nBlock,plotThings);

end