This matlab program can be used to analyze time traces of magnetic beads from a magnetic tweezer, and returns values for the bead radius and the magnetic force.

First setup the program by opening config.m. Here you insert values for the the requested configuration variables.

-zOffsetDataFile:
Insert the path to the trace containing a measurement of the z-direction with the magnet switched of or the magnet height sufficiently high or the magnetic field to be unnoticeable. This file should consist of three columns (X, Y, Z) in micrometers. X and Y values will not be used.

-zOffsetOutputFile:
Insert the path to save the file containing z-offset data to. If z-offsets were already saved insert here the path to the file containing the already saved files. This consists of one column containing the offsets in micrmeters.

-tracesFile:
Insert the path to the file containing the position trace of the bead. This file can include measurments at different magnet heights, as long as these heights are correctly represented in the magnetMotorFile (see below).
This file should consist of three columns (X, Y, Z) in micrometers. X and Y are perpendicular to the magnet force. Here X is in the short pendulum direction, and Y is in the long pendulum direction. Z is parallel to the magnet force.

-magnetMotorFile:
Insert the path to the file containing data on the magnet height, in one column. Units do not matter.

-sampleFreq:
Insert the acquisition frequency of your experiment.

-kT:
Insert the boltzmann constant times the temperature, in pN*nm.

-viscosity:
Insert the viscosity of the medium, in pN*s/nm^2.

-beadRadius:
Insert the radius of the used bead, as given by the manufacturer, in nm. The program will use this as an initial guess to fit the radius (and the force).

-plotThings:
true: plots are shown.
false: no plots are shown.

-nBlock:
The program divides the position traces in nBlock parts and averages over the nBlock Power Spectra. The largest possible number of blocks depends on the lengths of the individual plateaus. More blocks means faster and more accurate fitting. nBlock is a positive integer.

-zOffsetAlreadySaved:
true: skip the calculation of the z-direction offsets, because they were already calculated in an earlier run.
false: calculate the z-offsets (again).

-zOffsetAlreadySubtracted:
true: skip the calculation of the z-direction offsets and skip the subtraction of z-offsets, because in the data file they are already subtracted. (No valid paths will be needed for zOffsetDataFile and zOffsetOutputFile.)
false: calculate and subtract the z-offsets.

To run the program run main.m
