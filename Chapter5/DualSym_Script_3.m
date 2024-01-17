%% Script to simulate all the results of the first shape optimization problem
%% Add all folders
clear all
close all
clc
folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
rmpath(genpath('PrecomputedResults'));
% if ==0, then precomputed files will be used. Otherwise the simulation
% runs from scratch, what might take a long time.
perform_simulation = 0; 
%%
if perform_simulation == 1
    DielectricStraights3;
    PlotIterates3;
    PlotFinals3;
    Plot_Steps_vs_Measures3;
    Plot_Eps_vs_Measures3;
    Plot_k_vs_Measures3;
else
    addpath(genpath('PrecomputedResults')); % It will use the precomputed files now.
    PlotIterates3;
    PlotFinals3;
    Plot_Steps_vs_Measures3;
    Plot_EpsMu_vs_Measures3;
    Plot_k_vs_Measures3;
end

% rmpath(genpath(folder));