function [ZhongboFitShort, ZhongboForceShort] = ...
    analyzeZhongboShort(sampleFreq,extensionDNA,short,kT,beadRadius,viscosity,nBlock,plotThings)
%%% Analyzes the Power Spectral Density of a single trace in the short
%%% pendulum direction using Zhongbo Yu's method.

%%% Input: (sampleFreq,extensionDNA,short,beadRadius,kT,viscosity,nBlock,plotThings)
%%% - sampling frequency in Hz
%%% - extension of the DNA in nm
%%% - trace in short pendulum direction in nm
%%% - kT in pN nm
%%% - bead radius in nm
%%% - viscosity in pN s/nm^2
%%% - number of blocks used for blocked powerspectrum
%%% - show plots

%%% Output: [ZhongboFitShort, ZhongboForceShort]
%%% - results of the fit, force
%%

%%% Get an estimate of the force and the corner frequency
Fest = kT*extensionDNA/std(short)^2; %pN
cornerFreq = calcFcorner(Fest,extensionDNA,beadRadius,viscosity); %Hz
fitgood = cornerFreq < sampleFreq/2;


%%% Add a line to x to make the length 2^integer
short(end+1) = short(1);

%%% Calc PSD, find data points below 1/20 of fc (the corner frequency)
%%% Also throw away first point, which is merely mean(x)
%%% NOTE: 1/20 is hard-coded, seems to work fine for all data
[f, PSD, ~] = calcPSDBlock(short,sampleFreq,nBlock);
f(1) = []; PSD(1) = [];
goodinds = f > cornerFreq(1)/20;

%%% Max likelihood fit to PSD (using goodinds)
%%% According to the model by Lansdorp and Saleh (RSI 2012)
alpha0 = 1E-5; kappa0 = 4E-4;
[par, ~, ~] = fminsearch(@(par) costfunctionPSDfit(f(goodinds), par(1), par(2),sampleFreq, PSD(goodinds),kT), [alpha0 kappa0]);
alpha = par(1); kappa = abs(par(2)); ZhongboFitShort = [alpha kappa fitgood];
PSDmodel = analyticalPSDZhongboShort(alpha,kappa,sampleFreq,f,kT);
PSDmodelInit = analyticalPSDZhongboShort(alpha0,kappa0,sampleFreq,f,kT);
ZhongboForceShort = kappa*extensionDNA;

%%% Plot if plot
if plotThings;
    figure(4);
    loglog(f(goodinds),PSD(goodinds),'r-');
    hold on
    loglog(f(goodinds),PSDmodel(goodinds),'b-');
    loglog(f(goodinds),PSDmodelInit(goodinds),'y-');

    title('Fitting of Power Spectrum in short pendulum direction using Zhongbo Yus method');
    xlabel('frequency (Hz)');
    ylabel('Power Spectrum (nm^2/Hz)')
    legend(['Power spectrum, nBlock = ' num2str(nBlock)], 'Fitted analytical function','Initial guess for fit');
    hold off
end


end
