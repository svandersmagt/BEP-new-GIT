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

    %%% Write here some descriptive text about your data
    %%% ---
    %%% Example of force extension data
    %%% 21 kbp DNA, M270 beads, 1 mm gap verticaly oriented magnets
    tracesFile = configVariable.tracesFile;
    motorsFile = configVariable.magnetMotorFile;
    zOffsetsFile = configVariable.zOffsetOutputFile;

    %%% Read in data
    data = load(tracesFile);
    zmag = load(motorsFile);
    clear bead

    %%% For data = (long pendulum, short pendulum, z)
    if configVariable.pendulumOrder == 1;
        for i=1;
            bead(i).time = 1:length(data(:,1));
            bead(i).long = data(:,1)*1000; %nm
            bead(i).short = data(:,2)*1000; 
            bead(i).z = data(:,3)*1000;
        end
    end

    %%% For data = (short pendulum, long pendulum, z)
    if configVariable.pendulumOrder == 0;
        for i=1;
            bead(i).time = 1:length(data(:,1));
            bead(i).long = data(:,2)*1000; %nm
            bead(i).short = data(:,1)*1000; 
            bead(i).z = data(:,3)*1000;
        end
    end

    %%% Skip if offset is already subtracted
    if (offsetZ == 0);
        %%% Read in the previously determined z-offsets from file
        zOffData = load(zOffsetsFile);
        zOffsets = zOffData(:,1)*1000; %nm

        %%% Subtract previously determined z-offsets
        if length(zOffsets) == length(bead);
            for i=1;
                bead(i).z = bead(i).z - zOffsets(i);
            end
            display('Subtracted pre-determined z-offsets.')
        else
            display(['Number of beads in the big data file (' num2str(length(bead)) ...
                ') does not agree with the number of beads in the zoffset file ('...
                num2str(length(zOffsets)) ')'])
        end
    else
        display('Z-offsets were already subtracted')
    end
end