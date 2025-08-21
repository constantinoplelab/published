function PlotExtendedData8(datadir, codedir)
%PlotExtendedData8 - Plots Extended Data 4. 
% INPUTS:
%   datadir - Local directory where ProcessData_ExtendedData8.mat was saved after running ProcessData_ExtendedData8
%   codedir - Local directory where code (e.g., published/EstrousRPE/) was saved

%% Set up paths
s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codedir, 'IgnoreCase', ispc);
if ~onPath % only add code dir to path if it isn't already
    addpath(genpath(codedir))
end

%% Load data
load([datadir, 'ProcessData_ExtendedData8'], 'percentages')

%% PLOT %%
figure;
set(gcf, 'Color', [1 1 1], 'Units', 'Inches', 'Position', [0, 0, 17, 9],...
    'PaperUnits', 'Inches', 'PaperSize', [9, 5], 'renderer','painters')

%--------------------------------------------------------------------------
% ED8c. TH and ChR2 overlap
%--------------------------------------------------------------------
nexttile
histogram(percentages, binwidth=0.1, facecolor='k')
xlim([0 1.2])
xticks(0:0.1:1.2)
grid off
set(gca, 'TickDir', 'out'); box off;
xline(median(percentages), '-k'),
xline(1, '--k'),
ylabel('# Images')
xlabel('Overlap Area/Green Area')
title('c')
axis square

end



