clear all
close all
% clc
folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
rmpath(genpath('PrecomputedResults'));
material = 'silver';
tracking = 1;
Freq = 400 :100: 700;
DegreeOfVecSphHar = [2 4 6 8];

for jInd = 1 : 4
    given_frequency =Freq(jInd);
    OptiStraightFunction(given_frequency,DegreeOfVecSphHar, material,tracking)
end
%% Generate the silver table
fileID = fopen(strcat('SilverTable','.tex'),'w'); %open file for writing
fprintf(fileID, "\\begin{tabular}{llllll} \n");
fprintf(fileID, '\\hline\n');
fprintf(fileID, '&  \\multicolumn{4}{c}{\\textbf{Silver}} & \\\\ \\hline \n');
fprintf(fileID, '& $f [\\mathrm{THz}]$ & $400$ & $500$ & $600$& $700$ \\\\ \n');
fprintf(fileID, '& $\\lambda [\\mathrm{nm}]$ & $750$ & $600$ & $500$& $428$ \\\\ \\hline \\hline \n');

lenstrin = {'$\\frac{1}{4}\\lambda$', '$\\frac{1}{2}\\lambda$', '$\\lambda$', '$2\\lambda$'};
for klen = 1 : 4
    len = lenstrin{klen};
    for jInd = 1 : 2
%         jj
        if jInd == 1
            fprintf(fileID, strcat('{\\multirow{2}{*}{',len,'}} & $J_2$'));
            for kfreq = 1 : 4
                given_frequency =Freq(kfreq);
                load(strcat('Straights/HeightNo',num2str(klen),'_'...
                    ,num2str(given_frequency),'.mat'))
                pause(.1)
                fprintf(fileID, strcat('&',num2str(round(chir_alpha,2),'%10.2f')));
            end
            fprintf(fileID,'\\\\ \n');
        end
        if jInd == 2
            jInd
            fprintf(fileID, '& ﻿$J_{\\mathrm{HS}}$ ');
            for kfreq = 1 : 4
                given_frequency =Freq(kfreq);
                load(strcat('Straights/HeightNo',num2str(klen),'_'...
                    ,num2str(given_frequency),'.mat'))
                pause(.1)
                fprintf(fileID, strcat('&',num2str(round(smooth_relax_alpha,2),'%10.2f')));
            end
            fprintf(fileID,'\\\\ \\hline \n');
        end
    end
end
fprintf(fileID,'\\end{tabular}');
fclose(fileID);
%%
material = 'gold';
tracking = 1;
Freq = 400 :100: 700;
DegreeOfVecSphHar = [2 4 6 8];

for jInd = 1 : 4
    given_frequency =Freq(jInd);
    OptiStraightFunction(given_frequency,DegreeOfVecSphHar, material,tracking)
end
%%
fileID = fopen(strcat('GoldTable','.tex'),'w'); %open file for writing
fprintf(fileID, "\\begin{tabular}{llllll} \n");
fprintf(fileID, '\\hline\n');
fprintf(fileID, '&  \\multicolumn{4}{c}{\\textbf{Gold}} & \\\\ \\hline \n');
fprintf(fileID, '& $f [\\mathrm{THz}]$ & $400$ & $500$ & $600$& $700$ \\\\ \n');
fprintf(fileID, '& $\\lambda [\\mathrm{nm}]$ & $750$ & $600$ & $500$& $428$ \\\\ \\hline \\hline \n');

lenstrin = {'$\\frac{1}{4}\\lambda$', '$\\frac{1}{2}\\lambda$', '$\\lambda$', '$2\\lambda$'};
for klen = 1 : 4
    len = lenstrin{klen};
    for jInd = 1 : 2
%         jj
        if jInd == 1
            fprintf(fileID, strcat('{\\multirow{2}{*}{',len,'}} & $J_2$'));
            for kfreq = 1 : 4
                given_frequency =Freq(kfreq);
                load(strcat('Straights/HeightNo',num2str(klen),'_'...
                    ,num2str(given_frequency),'G.mat'))
                pause(.1)
                fprintf(fileID, strcat('&',num2str(round(chir_alpha,2),'%10.2f')));
            end
            fprintf(fileID,'\\\\ \n');
        end
        if jInd == 2
            jInd
            fprintf(fileID, '& ﻿$J_{\\mathrm{HS}}$ ');
            for kfreq = 1 : 4
                given_frequency =Freq(kfreq);
                load(strcat('Straights/HeightNo',num2str(klen),'_'...
                    ,num2str(given_frequency),'G.mat'))
                pause(.1)
                fprintf(fileID, strcat('&',num2str(round(smooth_relax_alpha,2),'%10.2f')));
            end
            fprintf(fileID,'\\\\ \\hline \n');
        end
    end
end
fprintf(fileID,'\\end{tabular}');
fclose(fileID);


