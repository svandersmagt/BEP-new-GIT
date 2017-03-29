
function fc = calcFcorner(forces,ext,R,eta)
% Runs a simple calculation to estimate the corner frequency in MT
% input
% forces in pN
% mean_z, extension in nm
% R, bead radius in nm
% eta, water viscosity in pN s/nm^2

fc = 1/(2*pi).* forces ./ ext ./(6*pi*eta*R);

end