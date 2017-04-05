function PSDmodel = analyticalPSDZhongboShort(alpha,kappa,sampleFreq,frequency,kT)
%%% Calculates the theoretical Power Spectral Density for Brownian motion
%%% in a harmonic well (Lansdorp and Saleh, RSI 2012)

%%% Input: (alpha,kappa,sampleFreq,frequency,kT)
%%% - alpha in pN s/nm
%%% - kappa in pN/nm
%%% - sampling frequency in Hz
%%% - frequency in Hz (vector)
%%% - kT in pN nm

%%% Output: PSDmodel in nm^2/Hz
%%
    preFactor = 4*kT*alpha/kappa^2;

    numerator = 2*alpha/kappa*sampleFreq*(sin(pi*frequency/sampleFreq).^2)*sinh(kappa/(alpha*sampleFreq));
    denominator = cos(2*pi*frequency/sampleFreq)-cosh(kappa/(alpha*sampleFreq));

    PSDmodel = preFactor*(1+numerator./denominator);
end