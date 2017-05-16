function [bead, zmag] = loadDataSubtractOffset(configVariable)
%%% STEP 2: LOAD DATA AND SUBTRACT Z-OFFSET
%%% -------------------------------------------------------------------
%%% If the z-offset is already subtracted in your data this can be
%%% indicated in the configuration file.

%%% Input: (configVariable)

%%% Output: [bead, zmag]
%%% - struct containing the time trace and position traces of the beads
%%% - vector containing motor data of the magnet height
%%
    format compact;
    offsetZ = configVariable.zOffsetAlreadySubtracted;
    skipColumn = configVariable.firstColumnIsTime;

    %%% Write here some descriptive text about your data
    %%% ---
    %%% Example of force extension data
    %%% 21 kbp DNA, M270 beads, 1 mm gap verticaly oriented magnets
    tracesFile = configVariable.tracesFile;
    if configVariable.motorDataInSameFile==0;
        motorsFile = configVariable.magnetMotorFile;
    end
    zOffsetsFile = configVariable.zOffsetOutputFile;

    %%% Read in data
    data = load(tracesFile);
    if configVariable.motorDataInSameFile==0;
        zmag = load(motorsFile);
    else
        zmag = data(:,skipColumn + 4);
    end
    clear bead

    %%% For data = (long pendulum, short pendulum, z)
    if configVariable.pendulumOrder == 1;
        bead.time = 1:length(data(:,1));
        bead.long = data(:,skipColumn + 1)*1000; %nm
        bead.short = data(:,skipColumn + 2)*1000; 
        bead.z = data(:,skipColumn + 3)*1000;
    end

    %%% For data = (short pendulum, long pendulum, z)
    if configVariable.pendulumOrder == 0;
        bead.time = 1:length(data(:,1));
        bead.long = data(:,skipColumn + 2)*1000; %nm
        bead.short = data(:,skipColumn + 1)*1000; 
        bead.z = data(:,skipColumn + 3)*1000;
    end

    %%% Skip if offset is already subtracted
    if (offsetZ == 0 && configVariable.zOffsetHighestPlateau == 0);
        %%% Read in the previously determined z-offsets from file
        zOffData = load(zOffsetsFile);
        zOffsets = zOffData(:,1)*1000; %nm

        %%% Subtract previously determined z-offsets
        if length(zOffsets) == length(bead);
            bead.z = bead.z - zOffsets;
            display('Subtracted pre-determined z-offsets.')
        else
            display(['Number of beads in the big data file (' num2str(length(bead)) ...
                ') does not agree with the number of beads in the zoffset file ('...
                num2str(length(zOffsets)) ')'])
        end
    elseif (offsetZ == 1 && configVariable.zOffsetHighestPlateau == 0)
        display('Z-offsets were already subtracted')
    elseif (configVariable.zOffsetHighestPlateau == 1)
        display('Z-offsets will be subtracted by using highest plateau')
    end
end