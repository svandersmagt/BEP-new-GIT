function [frequency,PowerSpectrum,maxTime] = calcPSDBlock(X,sampleFreq,nBlock)
%%% Calculates the Power Spectrum of a trace at a given sampling frequency
%%% by calculating the Power Spectrum of a given amount of blocks and averaging.

%%% Input: (X,sampleFreq,nBlock)
%%% - trace
%%% - sampling frequency in Hz
%%% - number of blocks to average over

%%% Output: [frequency,PowerSpectrum,maxTime]
%%% - vector of frequencies for Power Spectrum in Hz
%%% - Power Spectrum
%%% - duration of a block in s
%%
    lengthBlock = floor(length(X)/nBlock);
    Hannwindow = hann(lengthBlock)*sqrt(8/3);

    fNyq = sampleFreq/2;
    deltaT = 1/sampleFreq;
    time = (0 : deltaT : (lengthBlock-1)*deltaT); % same for every block
    maxTime = max(time,[],2);
    frequency = ((1 : lengthBlock) / maxTime)';
    ind = find(frequency <= fNyq); % only to the Nyquist f
    frequency = frequency(ind); % same for every block
    PSD = zeros(round(lengthBlock/2 -1),1);

    for i = 1:nBlock;
        FT = deltaT*fft(X((lengthBlock*(i-1)+1):lengthBlock*i).*Hannwindow);
        Ptemp = FT.*conj(FT)/maxTime;
        Ptemp = 2*Ptemp(ind);
        PSD = PSD(1:length(Ptemp))+Ptemp;

    end

    PowerSpectrum = PSD(1:length(Ptemp))/nBlock;
end
