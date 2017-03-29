function costfunc = costfunctionPSDfit(f,alpha,kappa,fs,PSD,kT)
% Costfunction for PSD Max Likelihood fit according to
% the model by Lansdorp and Saleh (RSI 2012)

eta = 1; % As far as I know, this is constant anyway, so why calc?
PSDmodel = analyticalPSDZhongboY(alpha,kappa,fs,f,kT);

costfunc = eta .* (PSD./PSDmodel + log(PSDmodel));
costfunc = sum(costfunc);

end