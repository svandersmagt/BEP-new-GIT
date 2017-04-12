%%% Force Analysis
%%% please read ReadMe and edit config.m

clear all; clc; close all; format compact;
configVariable = config();

%%% Step one: z-offset calculation
if(configVariable.zOffsetAlreadySaved == 0 && configVariable.zOffsetAlreadySubtracted == 0 && configVariable.zOffsetHighestPlateau == 0)
    display('Calculating z-offset')
    zOffsetCalculation(configVariable);
end

%%% Step two: loading data and subtracting offset
display('Loading bead data')
[bead, zmag] = loadDataSubtractOffset(configVariable);


%%% Step three: finding plateaus in magnet height
display('Finding plateaus')
[plat, zmags, nPlat] =  plateauFinding(zmag, configVariable);

%% start here to skip data loading

%%% This step only occurs if z-offsets were not subtracted yet
%%% Subtract the offset found by analyzing the highest plateau (with lowest
%%% force)
if configVariable.zOffsetHighestPlateau
   [bead] = subtractZOffsetHighestPlateau(zmags, plat, bead, configVariable); 
end

%%% Step four: analyze using different algorithms
display('Analyzing bead data')
[bead, forcesExponentialFit] = analyzeData(bead, nPlat, plat, zmags, configVariable);