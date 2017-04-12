function PSDmodel = analyticalPSDDaldropLong(force,sampleFreq,frequency,extensionDNA,L,beadRadius,kT,viscosity)
%%% Calculates the theoretical Power Spectral Density for Brownian motion
%%% in the long pendulum direction (Daldrop 2015)

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
    Cpar = (1-9/16*(1+L/beadRadius)^(-1)+1/8*(1+L/beadRadius)^(-3)...
        -45/256*(1+L/beadRadius)^(-4)-1/16*(1+L/beadRadius)^(-5))^(-1); %Daldrop eq(S10)
    Crot = 1 + 5/16*(1+L/beadRadius)^(-3); %Daldrop eq(S12)
    alphaX = 6*pi*viscosity*beadRadius*Cpar; %Daldrop eq(11)
    alphaPhi = 8*pi*viscosity*beadRadius^3*Crot; %Daldrop eq(11)

    fPlus = (force/extensionDNA*((extensionDNA+beadRadius)*beadRadius/(2*alphaPhi)...
        + 1/(2*alphaX) + 1/2*(((extensionDNA+beadRadius)*beadRadius/alphaPhi + 1/alphaX)^2 ...
        -4*extensionDNA*beadRadius/(alphaX*alphaPhi))^(1/2)))/(2*pi);
    
    fMin = (force/extensionDNA*((extensionDNA+beadRadius)*beadRadius/(2*alphaPhi)...
        + 1/(2*alphaX) - 1/2*(((extensionDNA+beadRadius)*beadRadius/alphaPhi + 1/alphaX)^2 ...
        -4*extensionDNA*beadRadius/(alphaX*alphaPhi))^(1/2)))/(2*pi); %f+ and f-, Daldrop eq(16)
    
    C = 2*pi*fPlus*extensionDNA/force - (extensionDNA+beadRadius)*beadRadius/alphaPhi; %a constant used in Daldrop eq(15)

    PSD1 = 4*kT/((2*pi)^2*(1+C^2*alphaX*alphaPhi/(beadRadius^2))); %first line of the analytical formula (eq(S13))
    PSD4 = zeros(length(frequency),1);
    
    for n=[-1 0];
        PSD2 = alphaPhi*C^2./(beadRadius^2)*1./(fPlus^2+abs(frequency+n*sampleFreq).^2) + 1./alphaX*1./(fMin^2+abs(frequency+n*sampleFreq).^2); %second line
        PSD3 = ((sin(pi/sampleFreq*abs(frequency+n*sampleFreq))).^2)./(pi/sampleFreq*abs(frequency+n*sampleFreq)).^2;
        PSD4 = PSD4 + PSD2.*PSD3;
    end



    PSDmodel = (PSD1*PSD4); %Daldrop eq (S13)

end