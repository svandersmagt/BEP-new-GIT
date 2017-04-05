function costFunction = costfunctionPSDfit(frequency,alpha,kappa,sampleFreq,PSD,kT)
% Costfunction for PSD Max Likelihood fit according to
% the model by Lansdorp and Saleh (RSI 2012)

eta = 1;
PSDmodel = analyticalPSDZhongboShort(alpha,kappa,sampleFreq,frequency,kT);

costFunction = eta .* (PSD./PSDmodel + log(PSDmodel));
costFunction = sum(costFunction);

end