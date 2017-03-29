function [DaldropXFit, DaldropXForce, DaldropXRadius] = analyzeDaldropX(fs,ext,x,y,R,kT,eta,nBlock,plot)
    %%% Analyzes the Power Spectral Density of a single trace in the x-direction,
    %%% using Daldrop's method.

    %%% Input: (fs,ext,x,y,R,kT,eta,nBlock,plot)
    %%% - fs in Hz, sampling frequency
    %%% - ext in nm, extension of the DNA (mean of z)
    %%% - x in nm, data trace in x direction, mean already subtracted
    %%% - y in nm, data trace in y direction, mean already subtracted
    %%% - R in nm, initial guess for the radius
    %%% - kT in pN nm
    %%% - eta, viscosity in pN s/nm^2
    %%% - nBlock, number of blocks used for blocked powerspectrum

    %%% Output: [DaldropXFit, DaldropXForce, DaldropXRadius]
    %%% - results of the fit

    %%

    %%% Get an estimate of the force and the corner frequency, using the y
    %%% direction, because x gives wrong values because it does not take
    %%% rotation into account
    Fest = kT*ext/std(y-mean(y))^2; %pN
    fc = calcFcorner(Fest,ext,R,eta); %Hz
    fitgood = fc < fs/2;

    %%% Add a line to x to make the length 2^integer
    x(end+1) = x(1);

    %%% Calc PSD, find data points below 1/20 of f_c (the corner frequency)
    %%% Also throw away first point, which is merely mean(x)
    %%% NOTE: 1/20 is hard-coded, seems to work fine for all data
    [f, PSD, ~] = calcPSDBlock(x,fs,nBlock);
    f(1) = []; PSD(1) = [];
    goodinds = f > fc(1)/20;
    f = f(goodinds);
    PSD = PSD(goodinds);
    
    %%% Weighted least squares fit to PSD
    %%% According to the model by Daldrop
    weight = (PSD * sqrt(nBlock));
    analyticalFunctionInverse = @(force, radius)(1./analyticalPSDDaldropX(force,fs,f,ext,radius,kT,eta));
    amplitudeBias = nBlock/(nBlock + 1);
    PSDInverse = 1./PSD;

    [par] = lsqnonlin(@(par) weight.*(PSDInverse - amplitudeBias*analyticalFunctionInverse(par(1),par(2))), [Fest,R]);
    DaldropXForce = par(1); DaldropXRadius = par(2); DaldropXFit = [DaldropXForce DaldropXRadius fitgood];
    
    PSDmodel = analyticalPSDDaldropX(DaldropXForce,fs,f,ext,DaldropXRadius,kT,eta);
    PSDmodelInit = analyticalPSDDaldropX(Fest,fs,f,ext,R,kT,eta);

    if plot;
        figure(5);
        loglog(f,PSD,'r-');
        hold on
        loglog(f,PSDmodel,'b-');
        loglog(f,PSDmodelInit,'y-');

        title('Fitting of Power Spectrum in x-direction using Daldrops method');
        xlabel('frequency (Hz)');
        ylabel('Power Spectrum (nm^2/Hz)')
        legend(['Power spectrum, nBlock = ' num2str(nBlock)], 'Fitted analytical function','Initial guess for fit');
        hold off
    end

end