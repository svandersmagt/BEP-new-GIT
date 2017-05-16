function [DaldropFitLong, DaldropForceLong, DaldropRadiusLong, cornerFreq, signal, errorLong] = ...
    analyzeDaldropLong(sampleFreq,extensionDNA,L,long,beadRadius,kT,viscosity,plotThings)
%%% Analyzes the Power Spectral Density of a single trace in the x-direction,
%%% using Daldrop's method.

%%% Input: (sampleFreq,extensionDNA,L,long,beadRadius,kT,viscosity,plotThings)
%%% - sampling frequency in Hz
%%% - extension of the DNA in nm
%%% - bead distance to surface
%%% - trace in long pendulum direction in nm
%%% - initial guess for the bead radius in nm
%%% - kT in pN nm
%%% - viscosity in pN s/nm^2
%%% - show plots

%%% Output: [DaldropFitLong, DaldropForceLong, DaldropRadiusLong, cornerFreq, signal, errorLong]
%%% - results of the fit, force and radius cornerfrequency and error
%%
    %%% Get an estimate of the force and the corner frequency, using the y
    %%% direction, because x gives wrong values because it does not take
    %%% rotation into account
    Fest = kT*(extensionDNA+beadRadius)/std(long)^2; %pN
    cornerFreq = calcFcorner(Fest,L,extensionDNA+beadRadius,beadRadius,viscosity); %Hz
    fitgood = cornerFreq < sampleFreq/2;
    
    %%% Add a line to x to make the length 2^integer
    long(end+1) = long(1);
    
    nBlock = round(length(long) / 5000) + 1;

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
    analyticalFunctionInverse = @(force, radius)(1./analyticalPSDDaldropLong(force,sampleFreq,f,extensionDNA,L,radius,kT,viscosity));
    amplitudeBias = nBlock/(nBlock + 1);
    PSDInverse = 1./PSD;

    [par,~,residual,~,~,~,jacobian] = lsqnonlin(@(par) weight.*(PSDInverse - amplitudeBias*analyticalFunctionInverse(par(1),par(2))), [Fest,beadRadius]);
    DaldropForceLong = par(1); DaldropRadiusLong = par(2); DaldropFitLong = fitgood;
    
    errorLong = nlparci(par',residual,'Jacobian',jacobian);
    
    PSDmodel = analyticalPSDDaldropLong(DaldropForceLong,sampleFreq,f,extensionDNA,L,DaldropRadiusLong,kT,viscosity);
    PSDmodelInit = analyticalPSDDaldropLong(Fest,sampleFreq,f,extensionDNA,L,beadRadius,kT,viscosity);
    signal = PSDmodel(1);

    if plotThings;
        figure(3);
        loglog(f,PSD,'r-');
        hold on
        loglog(f,PSDmodel,'b-');
        loglog(f,PSDmodelInit,'y-');
        loglog(f(1),signal,'kx');

        title('Fitting of Power Spectrum in long pendulum direction');
        xlabel('frequency (Hz)');
        ylabel('Power Spectrum (nm^2/Hz)')
        legend(['Power spectrum, nBlock = ' num2str(nBlock)], 'Fitted analytical function','Initial guess for fit');
        hold off
    end

end