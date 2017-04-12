This matlab program can be used to analyze time traces of magnetic beads from a magnetic tweezer, and returns values for the bead radius and the magnetic force.

--Setting up--
First setup the program by editing config.m. Insert values for the the requested configuration variables.

-zOffsetDataFile:
Insert the path to the trace containing a measurement of the z-direction with the magnet switched of or the magnet height sufficiently high or the magnetic field to be unnoticeable. This file should consist of three columns (X, Y, Z) in micrometers. X and Y values will not be used.

-zOffsetOutputFile:
Insert the path to save the file containing z-offset data to. If z-offsets were already saved insert here the path to the file containing the already saved files. This consists of one column containing the offsets in micrometers.

-tracesFile:
Insert the path to the file containing the position trace of the bead. This file can include measurments at different magnet heights, as long as these heights are correctly represented in the magnetMotorFile (see below).
This file should consist of three columns (X, Y, Z) in micrometers. X and Y are perpendicular to the magnet force. (See pendulumOrder) Z is parallel to the magnet force.

-magnetMotorFile:
Insert the path to the file containing data on the magnet height, in one column. Units do not matter. It is also possible to include the traces and motor information in a single file. See motorDataInSameFile. If this is the case magnetMotorFile does not have to be a valid path, because it will not be used.

-sampleFreq:
Insert the acquisition frequency of the experiment.

-pendulumOrder:
This variable is used to indicate the order of the different directions in the datafile. The input should be a boolean:
	-true: (X = long pendulum, Y = short pendulum, Z)
	-false: (X = short pendulum, Y = long pendulum, Z)

-firstColumnIsTime:
If the first column of the datafile is time, set to true. Else, set to false. Note: the time will still be computed with the given acquisition frequency. Setting to true will just skip the first column.

-motorDataInSameFile:
If the last column of the datafile is the motor data, set to true. Else, set to false.

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

-zOffsetHighestPlateau:
true: z-offset was not yet subtracted, but there is no separate offset file. The program will first find the magnet height plateaus and use the plateau with lowest force to calculate the z-offset.
false: z-offset subtraction will happen as indicated by the other variables.



--Running--
Running main.m starts the program. It consists of four main steps:

-Step 1: z-offset caclulation
Only runs if the z-offset was not yet calculated and a separate z-offset measurement was provided.
It saves the z-offset data to a file, indicated in config.m.

-Step 2: loading data and subtracting offset
Loads in the bead data and motor data from the indicated files. This can take a few minutes if the file is big. If necessary it subtracts the z-offset from the z-data.

-Step 3: finding plateaus in magnet height
Uses the magnet height motor data to find the plateaus where the magnet was at a constant height.

-Substep: use highest plateau to subtract z-offset
If indicated in config, the plateau with the lowest magnet force can be used to calculate and subtract the z-offset.

-Step 4: analyze the data
The powerspectrum of the position data is calculated, and in both the X and Y direction is fitted to force and radius. All results are saved to the struct 'bead'. Afterwards a double exponential is fitted to the force data to make a calibration curve for magnet height vs. force.



--Output--
-bead: struct containing data on the bead,
	-time: sample number
	-long: position trace, long pendulum direction
	-short: position trace, short pendulum direction
	-z: position trace, z-direction
	-extensionDNA: average DNA extension for every plateau
	-L: average z-height of the bead for every plateau
	-timeDomainForceLong/Short: force calculated in time domain direction, for every plateau
	-fitLong/Short: indicates for every plateau if the corner frequency is smaller than the Nyquist frequency. If not the fit is not good
	-forceLong/Short: fitted forces for every plateau
	-radiusLong/Short: fitted radii for every plateau
	-normalizedForce: Fshort/Flong for every plateau

-forcesExponentialFit: struct containing the double exponential fit parameters 	for both directions. These can be used to make a calibration curve:

	F(zMagnet) = delta + alpha0*exp(-zMagnet/zeta0) + alpha1*exp(-zMagnet/zeta1)

-nPlat: number of plateaus
-plat: struct containing the first and last index of every plateau
-zmag: magnet height data
-zmags: magnet height for every plateau



--Plots--
-Figure 1: trace used for finding z-offset, indicating lowest point
-Figure 2: magnet height vs time, indicating plateaus
-Figure 3: power spectrum, fit, and initial guess for fit in long pendulum direction. Shows all plateaus while running
-Figure 4: power spectrum, fit, and initial guess for fit in short pendulum direction. Shows all plateaus while running
-Figure 5: fitted force in the short pendulum direction relative to the fitted force in the long pendulum direction vs fitted force
-Figure 6: fitted radius in both directions vs fitted force
-Figure 7: magnet force vs magnet height, plus calibration curve fit