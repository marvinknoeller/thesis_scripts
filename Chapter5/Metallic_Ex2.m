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
    %%
    Example2
    Example2_Gold
    Example2_PlotStructures
    Example2_FrequencyScan
    Example2_FrequencyScan_Gold
else
    addpath(genpath('PrecomputedResults')); % It will use the precomputed files now.
    Example2_PlotStructures
    Example2_FrequencyScan
    Example2_FrequencyScan_Gold
end