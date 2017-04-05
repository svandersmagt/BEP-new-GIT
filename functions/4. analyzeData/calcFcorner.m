
function cornerFreq = calcFcorner(force,extensionDNA,beadRadius,viscosity)
%%% Runs a simple calculation to estimate the corner frequency in MT

%%% Input: (force, extensionDNA, beadRadius, viscosity)
%%% - magnet force in pN
%%% - extension of DNA in nm
%%% - estimate of the bead radius in nm
%%% - viscosity in pN s/nm^2

%%% Output: cornerFreq
%%% - corner frequency in Hz
%%
    cornerFreq = 1/(2*pi).* force ./ extensionDNA ./(6*pi*viscosity*beadRadius);
end