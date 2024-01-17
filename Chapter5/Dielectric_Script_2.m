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
    DielectricStraights2;
    PlotIterates2;
    PlotFinals2;
    Plot_Steps_vs_Measures2;
    Plot_Eps_vs_Measures2;
    Plot_k_vs_Measures2;
else
    addpath(genpath('PrecomputedResults')); % It will use the precomputed files now.
    PlotIterates2;
    PlotFinals2;
    Plot_Steps_vs_Measures2;
    Plot_eps_vs_Measures2;
    Plot_k_vs_Measures2;
end

% rmpath(genpath(folder));