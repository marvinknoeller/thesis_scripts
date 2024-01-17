clear all
clc
close all
folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
rmpath(genpath('PrecomputedResults'));
%%
perform_simulation = 0; 
%%
if perform_simulation == 1
    Example1
    Example1_PlotTubes
    Example1_FrequencyScan
    Example1_FrequencyScan_Gold
else
    addpath(genpath('PrecomputedResults')); % It will use the precomputed files now.
    Example1_PlotTubes
    Example1_FrequencyScan
    Example1_FrequencyScan_Gold
end
