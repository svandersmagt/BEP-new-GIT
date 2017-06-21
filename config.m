function [configVariable] = config()

%%% For calculation of the z-offset, data without a magnetic field is
%%% needed.
%%% File name:
configVariable.zOffsetDataFile = '..\BEPmaarnietgit\offset.txt';
%%% Z-offset results are saved to a file.
%%% File name:
configVariable.zOffsetOutputFile = '..\BEPmaarnietgit\FX_offsets.txt';


%%% Data file containing x, y, z traces of beads:
configVariable.tracesFile = '..\BEPmaarnietgit\WithoutReferenceCorrection\RdRp_4.txt';
%%% Data file containing motor data for the magnet:
configVariable.magnetMotorFile = '..\BEPmaarnietgit\bead_motors.txt';
%%% Acquisition frequency:
configVariable.sampleFreq = 114; %acquisition frequency in Hz

%%% Pendulum geometry
%%% true: data = (Long pendulum, Short pendulum, Z) 
%%% false: data = (Short pendulum, Long pendulum, Z)
configVariable.pendulumOrder = true;
%%% First column is time? (note: time will not be used as knowing the
%%% sampling frequency makes it obsolete, but the time column has to be
%%% skipped)
configVariable.firstColumnIsTime = true;
%%% Motor data in same file?
configVariable.motorDataInSameFile = true;

%%% Constants
configVariable.kT = 4.1; %pN nm
configVariable.viscosity = 10E-10; %viscosity in pN s/nm^2
configVariable.beadRadius = 1400; %bead radius in nm, as given by the manufacturer

%%% Options:
%%% Plot things:
configVariable.plotThings = true;
configVariable.plotSpectrum = false;

configVariable.maxMagnetHeigth = 20;

%%% Z-offsets are already saved to the file specified above:
configVariable.zOffsetAlreadySaved = false;

%%% There are no Z-offsets because they are already subtracted:
configVariable.zOffsetAlreadySubtracted = false;

%%% Subtract z-offsets using the highest plateau
configVariable.zOffsetHighestPlateau = true;
end