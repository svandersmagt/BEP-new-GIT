
function cornerFreq = calcFcorner(force,L,extensionDNA,beadRadius,viscosity)
%%% Runs a simple calculation to estimate the corner frequency in MT

%%% Input: (force, extensionDNA, beadRadius, viscosity)
%%% - magnet force in pN
%%% - extension of DNA in nm
%%% - estimate of the bead radius in nm
%%% - viscosity in pN s/nm^2

%%% Output: cornerFreq
%%% - corner frequency in Hz
%%
    %cornerFreq = 1/(2*pi).* force ./ extensionDNA ./(6*pi*viscosity*beadRadius);
    
    Cpar = (1-9/16*(1+L/beadRadius).^(-1)+1/8*(1+L/beadRadius).^(-3) ...
        -45/256*(1+L/beadRadius).^(-4)-1/16*(1+L/beadRadius).^(-5)).^(-1); %Daldrop eq(S10)
    alphaY = 6*pi*viscosity*beadRadius*Cpar; %Daldrop eq(6)
    kappa = force./extensionDNA;
    cornerFreq = kappa./(2*pi*alphaY);
end