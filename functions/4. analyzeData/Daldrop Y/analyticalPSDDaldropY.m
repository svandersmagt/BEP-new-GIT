function PSDmodel = analyticalPSDDaldropY(Fmag,fs,f,L,R,kT,eta)
%%% Calculates the theoretical Power Spectral Density for Brownian motion
%%% in the y-direction (Daldrop 2015)
%%% Input: Fmag in pN (number)
%%%        fs in Hz (number)
%%%        f in Hz (vector)
%%%        L in nm (number)
%%%        R in nm (number)
%%%        kT in pN nm
%%%        eta in pN s / nm

%%% Output: PSDmodel in nm^2/Hz (vector)
%%

Cpar = (1-9/16*(1+L/R)^(-1)+1/8*(1+L/R)^(-3)-45/256*(1+L/R)^(-4)-1/16*(1+L/R)^(-5))^(-1); %Daldrop eq(S10)
alphaY = 6*pi*eta*R*Cpar; %Daldrop eq(6)
kappa = Fmag/L;

fc = kappa/(2*pi*alphaY);

PSD1 = 4*kT*alphaY/(kappa)^2; %Daldrop eq(5)
PSD4 = zeros(length(f),1);
for n=[-1 0];
    PSD2 = 1./(1+(abs(f+n*fs)./fc).^2);
    PSD3 = ((sin(pi/fs*abs(f+n*fs))).^2)./(pi/fs*abs(f+n*fs)).^2;
    PSD4 = PSD4 + PSD2.*PSD3;
end



PSDmodel = PSD1*PSD4; %Daldrop eq (5)

end