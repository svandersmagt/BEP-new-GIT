function [ZhongboYFit, ZhongboYForce] = analyzeZhongboY(fs,ext,y,kT,R,eta,nBlock,plot)
%%% Analyzes the Power Spectral Density of a single trace in the y-direction,
%%% using Zhongbo's method.

%%% Input: (fs,ext,y,kT,R,eta,nBlock,plot)
%%% - fs in Hz, sampling frequency
%%% - ext in nm, extension of the DNA (mean of z)
%%% - y in nm, data trace in y direction, mean already subtracted
%%% - kT in pN nm
%%% - R in nm, given radius of bead (in this method the radius is not fitted)
%%% - eta, viscosity in pN s/nm^2
%%% - nBlock, number of blocks used for blocked powerspectrum

%%% Output: [ZhongboYFit, ZhongboYForce]
%%% - results of the fit

%%

%%% Get an estimate of the force and the corner frequency
Fest = kT*ext/std(y)^2; %pN
fc = calcFcorner(Fest,ext,R,eta); %Hz
fitgood = fc < fs/2;


%%% Add a line to x to make the length 2^integer
y(end+1) = y(1);

%%% Calc PSD, find data points below 1/20 of fc (the corner frequency)
%%% Also throw away first point, which is merely mean(x)
%%% NOTE: 1/20 is hard-coded, seems to work fine for all data
[f, PSD, ~] = calcPSDBlock(y,fs,nBlock);
f(1) = []; PSD(1) = [];
goodinds = f > fc(1)/20;

%%% Max likelihood fit to PSD (using goodinds)
%%% According to the model by Lansdorp and Saleh (RSI 2012)
alpha0 = 1E-5; kappa0 = 4E-4;
[par, ~, ~] = fminsearch(@(par) costfunctionPSDfit(f(goodinds), par(1), par(2),fs, PSD(goodinds),kT), [alpha0 kappa0]);
alpha = par(1); kappa = abs(par(2)); ZhongboYFit = [alpha kappa fitgood];
PSDmodel = analyticalPSDZhongboY(alpha,kappa,fs,f,kT);
PSDmodelInit = analyticalPSDZhongboY(alpha0,kappa0,fs,f,kT);
ZhongboYForce = kappa*ext;

%%% Plot if plot
if plot;
    figure(4);
    loglog(f(goodinds),PSD(goodinds),'r-');
    hold on
    loglog(f(goodinds),PSDmodel(goodinds),'b-');
    loglog(f(goodinds),PSDmodelInit(goodinds),'y-');

    title('Fitting of Power Spectrum in y-direction using Zhongbo Yus method');
    xlabel('frequency (Hz)');
    ylabel('Power Spectrum (nm^2/Hz)')
    legend(['Power spectrum, nBlock = ' num2str(nBlock)], 'Fitted analytical function','Initial guess for fit');
    hold off
end


end
