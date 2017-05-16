function [extensionDNA, L, timeDomainForceLong, timeDomainForceShort, ...
    DaldropFitLong, DaldropForceLong, DaldropRadiusLong, cornerFreqLong, signalLong, errorLong,...
    DaldropFitShort, DaldropForceShort, DaldropRadiusShort, cornerFreqShort, signalShort, errorShort]...
    = analyzeTrace(long, short, z, sampleFreq, beadRadius, kT, viscosity, plotThings)
%%% Function to analyze magnetic tweezers time traces in three different ways:
%%% -Daldrop's method of fitting the spectrum to the analytical solution in the short pendulum direction.
%%% -Daldrop's method of fitting the spectrum to the analytical solution in the long pendulum direction, taking bead rotations into account.

%%% Input: (long, short, z, sampleFreq, beadRadius, kT, viscosity, plotThings)
%%% - trace in long pendulum direction in nm
%%% - trace in short pendulum direction in nm
%%% - trace z in nm
%%% - sampling frequency in Hz
%%% - initial guess for the bead Radius in nm
%%% - kT in pN nm
%%% - viscosity in pN s/nm^2
%%% - show plots

%%% Output: [extensionDNA, L, timeDomainForceLong, timeDomainForceShort, ...
%%% DaldropFitLong, DaldropForceLong, DaldropRadiusLong, cornerFreqLong, signalLong, errorLong,...
%%% DaldropFitShort, DaldropForceShort, DaldropRadiusShort, cornerFreqShort, signalShort, errorShort]
%%% - results of analysis in time domain, and using Daldrop's methods
%%% including errors
%%
    %%% Analyze the real time fluctuations
    meanLong = mean(long);
    meanShort = mean(short);
    meanZ = mean(z);
    stdLong  = std(long);
    stdShort  = std(short); 
    
    if(stdLong^2<stdShort^2)
        display('Warning: the variance of the position in the long pendulum') 
        display('direction is smaller than that of the short direction. It is') 
        display('likely the two pendulum directions were switched. See config.m')
        pause;
    end
    
    extensionDNA = mean(sqrt((long-meanLong).^2 + (short-meanShort).^2 + z.^2));
    timeDomainForceLong = kT*(extensionDNA+beadRadius)./stdLong^2;
    timeDomainForceShort = kT*extensionDNA./stdShort^2;
    L = meanZ;
   

    %%% The three methods are implemented in their own functions:
    fprintf('Daldrop Long Pendulum Fit \n');
    [DaldropFitLong, DaldropForceLong, DaldropRadiusLong, cornerFreqLong, signalLong, errorLong] = analyzeDaldropLong(sampleFreq,...
        extensionDNA,L,long-meanLong,beadRadius,kT,viscosity,plotThings);
    
    fprintf('Daldrop Short Pendulum Fit \n');
    [DaldropFitShort, DaldropForceShort, DaldropRadiusShort, cornerFreqShort, signalShort, errorShort] = analyzeDaldropShort(...
        sampleFreq,extensionDNA,L,short-meanShort,beadRadius,kT,viscosity,plotThings);
end