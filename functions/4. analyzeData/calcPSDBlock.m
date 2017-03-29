function [f,P,T] = calcPSDBlock(X,sampling_f,nBlock)
% function [f,P,T] = calc_powersp(X,sampling_f)
% Calculates the powerspectrum P as function of frequency f
% of imput data X sampled at frequency sampling_f

lengthBlock = floor(length(X)/nBlock);
Hannwindow = hann(lengthBlock)*sqrt(8/3);

fNyq    =   sampling_f / 2;
delta_t =   1 / sampling_f;
time    =   [0 : delta_t : (lengthBlock-1)*delta_t]';
T       =   max(time);
f       =   ([1 : lengthBlock] / T)';
ind     =   find(f <= fNyq); % only to the Nyquist f
f       =   f(ind);
PSD     =   zeros(round(lengthBlock/2 -1),1);


for i = 1:nBlock;
    FT      =   delta_t*fft(X((lengthBlock*(i-1)+1):lengthBlock*i).*Hannwindow);
    Ptemp       =   FT .* conj(FT) / T;
    Ptemp       =   2 * Ptemp(ind);
    PSD     =   PSD(1:length(Ptemp)) + Ptemp;
    
end
P      =   PSD(1:length(Ptemp))/nBlock;
%P = Ptemp.*nBlock/(nBlock + 1);
end
