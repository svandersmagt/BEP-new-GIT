clear all; clc; close all; format compact;
configVariable = config();

%%% Step one: z-offset calculation
if(configVariable.zOffsetAlreadySaved == 0 && configVariable.zOffsetAlreadySubtracted == 0)
    display('Calculating z-offset')
    zOffsetCalculation(configVariable);
end

%%% Step two: loading data and subtracting offset
display('Loading bead data')
[bead, zmag] = loadDataSubtractOffset(configVariable);

%%% Step three: finding plateaus in magnet height
display('Finding plateaus')
[plat, zmags, nPlat] =  plateauFinding(zmag, configVariable);

%%% Step four: analyze using different algorithms
display('Analyzing bead data')
[bead, forcesExponentialFit] = analyzeData(bead, nPlat, plat, zmags, configVariable);