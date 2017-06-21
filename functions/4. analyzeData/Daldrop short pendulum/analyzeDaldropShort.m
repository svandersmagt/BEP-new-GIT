function [DaldropFitShort, DaldropForceShort, DaldropRadiusShort, cornerFreq, signal, errorShort] = ...
    analyzeDaldropShort(sampleFreq,extensionDNA,L,short,beadRadius,kT,viscosity,plotThings)
%%% Analyzes the Power Spectral Density of a single trace in the short
%%% pendulum direction using Daldrop's method.

%%% Input: (sampleFreq,extensionDNA,L,short,beadRadius,kT,viscosity,plotThings)
%%% - sampling frequency in Hz
%%% - extension of the DNA in nm
%%% - bead distance to surface in nm
%%% - trace in short pendulum direction in nm
%%% - initial guess for the bead radius in nm
%%% - kT in pN nm
%%% - viscosity in pN s/nm^2
%%% - show plots

%%% Output: [DaldropFitShort, DaldropForceShort, DaldropRadiusShort, cornerFreq, signal, errorShort]
%%% - results of the fit, force and radius cornerfrequency and error
%%
    %%% Get an estimate of the force and the corner frequency
    Fest = kT*extensionDNA/std(short-mean(short))^2; %pN
    cornerFreq = calcFcorner(Fest,L,extensionDNA,beadRadius,viscosity); %Hz
    fitgood = cornerFreq < sampleFreq/2;

    %%% Add a line to x to make the length 2^integer
    short(end+1) = short(1);
    
    nBlock = round(length(short) / 5000) + 1;

    %%% Calc PSD, find data points below 1/20 of f_c (the corner frequency)
    %%% Also throw away first point, which is merely mean(x)
    %%% NOTE: 1/20 is hard-coded, seems to work fine for all data
    [f, PSD, ~] = calcPSDBlock(short,sampleFreq,nBlock);
    f(1) = []; PSD(1) = [];
    goodinds = f > cornerFreq(1)/20;
    f = f(goodinds);
    PSD = PSD(goodinds);

    %%% Weighted least squares fit to PSD
    %%% According to the model by Daldrop
    weight = (PSD * sqrt(nBlock));
    analyticalFunctionInverse = @(force, radius)(1./analyticalPSDDaldropShort(force,sampleFreq,f,extensionDNA,L,radius,kT,viscosity));
    amplitudeBias = nBlock/(nBlock + 1);
    PSDInverse = 1./PSD;

    [par,~,residual,~,~,~,jacobian] = lsqnonlin(@(par) weight.*(PSDInverse - amplitudeBias*analyticalFunctionInverse(par(1),par(2))), [Fest,beadRadius]);
    DaldropForceShort = par(1); DaldropRadiusShort = par(2); DaldropFitShort = fitgood;
    
    errorShort = nlparci(par',residual,'Jacobian',jacobian);

    PSDmodel = analyticalPSDDaldropShort(DaldropForceShort,sampleFreq,f,extensionDNA,L,DaldropRadiusShort,kT,viscosity);
    PSDmodelInit = analyticalPSDDaldropShort(Fest,sampleFreq,f,extensionDNA,L,beadRadius,kT,viscosity);
    signal = PSDmodel(1);

    if plotThings;
        figure(4);
        loglog(f,PSD,'r-');
        hold on
        loglog(f,PSDmodel,'b-');
        loglog(f,PSDmodelInit,'y-');

        title('Fitting of Power Spectrum in short pendulum direction');
        xlabel('frequency (Hz)');
        ylabel('Power Spectrum (nm^2/Hz)');
        legend(['Power spectrum, nBlock = ' num2str(nBlock)], 'Fitted analytical function','Initial guess for fit');
        hold off
    end

end