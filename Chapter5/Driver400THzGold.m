
clear all
clc
close all
folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
%%
% This will run 100 simulations of the free form optimization at 400THz.
% Afterwards, a sample needs to be chosen manually. This is necessary since
% the optimization does not respect the fact that the object may not
% intersect itself.
% If the sample is chosen in the file Example3_PlotStructures and
% Example3_frequencyScan400Gold, then you obtain both the final structure and
% the frequency scan.
NExample3_400THz_Gold

%%
% Running this code without any edits will provide you with the results for
% the sample that I used in my thesis.
% Example3_frequencyScan400Gold
