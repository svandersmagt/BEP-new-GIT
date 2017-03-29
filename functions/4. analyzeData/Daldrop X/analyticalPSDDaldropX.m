function PSDmodel = analyticalPSDDaldropX(Fmag,fs,f,L,R,kT,eta)
%%% Calculates the theoretical Power Spectral Density for Brownian motion
%%% in the x-direction (Daldrop 2015)
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
Crot = 1 + 5/16*(1+L/R)^(-3); %Daldrop eq(S12)
alphaX = 6*pi*eta*R*Cpar;% + 8*pi*eta*R*Crot/(1+L/R)^2; %Daldrop eq(7) or Daldrop eq(11), it is not clear to me which one should be used.. eq(11) seems to give better results
alphaPhi = 8*pi*eta*R^3*Crot; %Daldrop eq(11)

fPlus = (Fmag/L*((L+R)*R/(2*alphaPhi) + 1/(2*alphaX) + 1/2*(((L+R)*R/alphaPhi + 1/alphaX)^2-4*L*R/(alphaX*alphaPhi))^(1/2)))/(2*pi);
fMin = (Fmag/L*((L+R)*R/(2*alphaPhi) + 1/(2*alphaX) - 1/2*(((L+R)*R/alphaPhi + 1/alphaX)^2-4*L*R/(alphaX*alphaPhi))^(1/2)))/(2*pi); %f+ and f-, Daldrop eq(16)
C = 2*pi*fPlus*L/Fmag - (L+R)*R/alphaPhi; %a constant used in Daldrop eq(15)

PSD1 = 4*kT/((2*pi)^2*(1+C^2*alphaX*alphaPhi/(R^2))); %first line of the analytical formula
PSD4 = zeros(length(f),1);
for n=[-1 0];
    PSD2 = alphaPhi*C^2./(R^2)*1./(fPlus^2+abs(f+n*fs).^2) + 1./alphaX*1./(fMin^2+abs(f+n*fs).^2); %second line
    PSD3 = ((sin(pi/fs*abs(f+n*fs))).^2)./(pi/fs*abs(f+n*fs)).^2;
    PSD4 = PSD4 + PSD2.*PSD3;
end



PSDmodel = (PSD1*PSD4); %Daldrop eq (S13)

end