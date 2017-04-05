function PSDmodel = analyticalPSDDaldropShort(force,sampleFreq,frequency,extensionDNA,beadRadius,kT,viscosity)
%%% Calculates the theoretical Power Spectral Density for Brownian motion
%%% in the short pendulum direction (Daldrop 2015)

%%% Input: (force,sampleFreq,frequency,extensionDNA,beadRadius,kT,viscosity) 
%%% - magnet force in pN
%%% - sampling frequency in Hz
%%% - frequencies for powerspectrum in Hz (vector)
%%% - DNA extension in nm
%%% - bead radius in nm
%%% - kT in pN nm
%%% - viscosity in pN s/nm^2

%%% Output: PSDmodel in nm^2/Hz
%%
    Cpar = (1-9/16*(1+extensionDNA/beadRadius)^(-1)+1/8*(1+extensionDNA/beadRadius)^(-3) ...
        -45/256*(1+extensionDNA/beadRadius)^(-4)-1/16*(1+extensionDNA/beadRadius)^(-5))^(-1); %Daldrop eq(S10)
    alphaY = 6*pi*viscosity*beadRadius*Cpar; %Daldrop eq(6)
    kappa = force/extensionDNA;

    cornerFreq = kappa/(2*pi*alphaY);

    PSD1 = 4*kT*alphaY/(kappa)^2; %Daldrop eq(5)
    PSD4 = zeros(length(frequency),1);
    for n=[-1 0];
        PSD2 = 1./(1+(abs(frequency+n*sampleFreq)./cornerFreq).^2);
        PSD3 = ((sin(pi/sampleFreq*abs(frequency+n*sampleFreq))).^2)./(pi/sampleFreq*abs(frequency+n*sampleFreq)).^2;
        PSD4 = PSD4 + PSD2.*PSD3;
    end

    PSDmodel = PSD1*PSD4; %Daldrop eq (5)
end