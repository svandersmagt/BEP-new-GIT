function [bead, zmag] = loadDataSubtractOffset(configVariable)
%%% STEP 2: LOAD DATA AND SUBTRACT Z-OFFSET

format compact;
offsetZ = configVariable.zOffsetAlreadySubtracted;

%%% Write here some descriptive text about your data
%%% ---
%%% Example of force extension data
%%% 21 kbp DNA, M270 beads, 1 mm gap verticaly oriented magnets
tracesFile = '..\BEPmaarnietgit\bead.txt';
motorsFile = '..\BEPmaarnietgit\bead_motors.txt';
zOffsetsFile = '..\BEPmaarnietgit\FX_offsets.txt';

%%% Read in data
data = load(tracesFile);
zmag = load(motorsFile);

for i=1;
    bead(i).time = 1:length(data(:,1));
    bead(i).x = data(:,1)*1000; %nm
    bead(i).y = data(:,1)*1000; %only one trace in this example, so x and y are the same
    bead(i).z = data(:,2)*1000;
end

if (offsetZ == 0);
    %%% Read in the previously determined z-offsets from file
    zoff_data = load(zOffsetsFile);
    zoffsets = zoff_data(:,2)*1000; %nm
    
    %%% Subtract previously determined z-offsets
    if length(zoffsets) == length(bead);
        for i=1;
            bead(i).z = bead(i).z - zoffsets(i);
        end
        display('Subtracted pre-determined z-offsets.')
    else
        display(['Number of beads in the big data file (' num2str(length(bead)) ') does not agree with the number of beads in the zoffset file (' num2str(length(zoffsets)) ')'])
    end
else
    display('Z-offsets were already subtracted')
end
end