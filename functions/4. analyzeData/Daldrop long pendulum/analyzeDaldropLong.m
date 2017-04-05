function [DaldropFitLong, DaldropForceLong, DaldropRadiusLong] = ...
    analyzeDaldropLong(sampleFreq,extensionDNA,long,short,beadRadius,kT,viscosity,nBlock,plotThings)
%%% Analyzes the Power Spectral Density of a single trace in the x-direction,
%%% using Daldrop's method.

%%% Input: (sampleFreq,extension,long,short,beadRadius,kT,viscosity,nBlock,plotThings)
%%% - sampling frequency in Hz
%%% - extension of the DNA in nm
%%% - trace in long pendulum direction in nm
%%% - trace in short pendulum direction in nm
%%% - initial guess for the bead radius in nm
%%% - kT in pN nm
%%% - viscosity in pN s/nm^2
%%% - number of blocks used for blocked powerspectrum
%%% - show plots

%%% Output: [DaldropFitLong, DaldropForceLong, DaldropRadiusLong]
%%% - results of the fit, force and radius
%%
    %%% Get an estimate of the force and the corner frequency, using the y
    %%% direction, because x gives wrong values because it does not take
    %%% rotation into account
    Fest = kT*extensionDNA/std(short-mean(short))^2; %pN
    cornerFreq = calcFcorner(Fest,extensionDNA,beadRadius,viscosity); %Hz
    fitgood = cornerFreq < sampleFreq/2;

    %%% Add a line to x to make the length 2^integer
    long(end+1) = long(1);

    %%% Calc PSD, find data points below 1/20 of f_c (the corner frequency)
    %%% Also throw away first point, which is merely mean(x)
    %%% NOTE: 1/20 is hard-coded, seems to work fine for all data
    [f, PSD, ~] = calcPSDBlock(long,sampleFreq,nBlock);
    f(1) = []; PSD(1) = [];
    goodinds = f > cornerFreq(1)/20;
    f = f(goodinds);
    PSD = PSD(goodinds);
    
    %%% Weighted least squares fit to PSD
    %%% According to the model by Daldrop
    weight = (PSD * sqrt(nBlock));
    analyticalFunctionInverse = @(force, radius)(1./analyticalPSDDaldropLong(force,sampleFreq,f,extensionDNA,radius,kT,viscosity));
    amplitudeBias = nBlock/(nBlock + 1);
    PSDInverse = 1./PSD;

    [par] = lsqnonlin(@(par) weight.*(PSDInverse - amplitudeBias*analyticalFunctionInverse(par(1),par(2))), [Fest,beadRadius]);
    DaldropForceLong = par(1); DaldropRadiusLong = par(2); DaldropFitLong = [DaldropForceLong DaldropRadiusLong fitgood];
    
    PSDmodel = analyticalPSDDaldropLong(DaldropForceLong,sampleFreq,f,extensionDNA,DaldropRadiusLong,kT,viscosity);
    PSDmodelInit = analyticalPSDDaldropLong(Fest,sampleFreq,f,extensionDNA,beadRadius,kT,viscosity);

    if plotThings;
        figure(5);
        loglog(f,PSD,'r-');
        hold on
        loglog(f,PSDmodel,'b-');
        loglog(f,PSDmodelInit,'y-');

        title('Fitting of Power Spectrum in long pendulum direction using Daldrops method');
        xlabel('frequency (Hz)');
        ylabel('Power Spectrum (nm^2/Hz)')
        legend(['Power spectrum, nBlock = ' num2str(nBlock)], 'Fitted analytical function','Initial guess for fit');
        hold off
    end

end