function [configVariable] = config()

%%% For calculation of the z-offset, data without a magnetic field is
%%% needed.
%%% File name:
configVariable.zOffsetDataFile = '..\BEPmaarnietgit\offset.txt';

%%% Z-offset results are saved to a file.
%%% File name:
configVariable.zOffsetOutputFile = '..\BEPmaarnietgit\FX_offsets.txt';


%%% Data file containing x, y, z traces of beads:
configVariable.tracesFile = '..\BEPmaarnietgit\DaldropData.txt';
%%% Data file containing motor data for the magnet:
configVariable.magnetMotorFile = '..\BEPmaarnietgit\DaldropMotors.txt';
%%% Acquisition frequency:
configVariable.sampleFreq = 2800; %acquisition frequency in Hz
%%% true: (Long pendulum, Short pendulum, Z) 
%%% false:(Short pendulum, Long pendulum, Z)
configVariable.pendulumOrder = true;

%%% Constants
configVariable.kT = 4.1; %pN nm
configVariable.viscosity = 10E-10; %viscosity in pN s/nm^2
configVariable.beadRadius = 515; %bead radius in nm, as given by the manufacturer

%%% Options:
%%% Plot things:
configVariable.plotThings = true;

%%% Number of blocks for blocked power spectrum calculation:
%%% Largest possible number of blocks depends on your trace length, more
%%% blocks speeds up the fitting, and might improve accuracy.
configVariable.nBlock = 40;

%%% Z-offsets are already saved to the file specified above:
configVariable.zOffsetAlreadySaved = false;

%%% There are no Z-offsets because they are already subtracted:
configVariable.zOffsetAlreadySubtracted = true;
end