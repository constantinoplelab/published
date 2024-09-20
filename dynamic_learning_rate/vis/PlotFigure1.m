function [] = PlotFigure1(datadir, codedir)
%PlotFigure1 - Plots Figure 1. Must be run after ProcessData_Figure1.
% INPUTS:
%   datadir - Directory where output from ProcessData_Figure1 was saved
%   codedir - Directory of code (e.g., published/dynamic_learning_rate/)

%% Set Path and Load Data

s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codedir, 'IgnoreCase', ispc);
if ~onPath % only add code dir to path if it already isn't
    addpath(genpath(codedir))
end

% Load pre-processed data. 
% To see pipeline for Figure 1, see Process/ProcessData_Figure1.m or 
% Process/ProcessPaperData.m
load([datadir, 'Figure1_Data'],...
    'ratList', 'exampleRat', 'exampleRatI', 'itiByBlk', 'pITIbyBlkRat',...
    'betasAll', 'pITIbyBlkPop', 'pRegressionPop');


%%
clear l
ht = 1;
wd = 4;

figure

%--------------------------------------------------------------------------
% Block
%--------------------------------------------------------------------------
subplot(ht, wd, 1)
title('Block')

%--------------------------------------------------------------------------
% Example rat ITI by blk
%--------------------------------------------------------------------------

subplot(ht, wd, 2); hold on
errorbar(1, itiByBlk.raw(exampleRatI, 1), itiByBlk.raw(exampleRatI, 2),...
    'bo', markerfacecolor='b', capsize=0, linewidth=2)
errorbar(2, itiByBlk.raw(exampleRatI, 3), itiByBlk.raw(exampleRatI, 4),...
    'ko', markerfacecolor='k', capsize=0, linewidth=2)
errorbar(3, itiByBlk.raw(exampleRatI, 5), itiByBlk.raw(exampleRatI, 6),...
    'ro', markerfacecolor='r', capsize=0, linewidth=2)
xlim([0.5 3.5])
ylim([2.7 5])

xticks(1:2)
xticklabels({'Low', 'Mixed', 'High'})
xlabel('Block')
ylabel('Trial initiation time (s)')
title({exampleRat, pITIbyBlkRat(exampleRatI,:)})

%--------------------------------------------------------------------------
% Population ITI by blk
%--------------------------------------------------------------------------

subplot(ht, wd, 3); hold on

plot(1:3, itiByBlk.z(:,[1, 3, 5]), color=[0.7 0.7 0.7 0.15])

errorbar(1,...
    mean(itiByBlk.z(:,1)),...
    std(itiByBlk.z(:,1)),...
    'bo', markerfacecolor='b', linewidth=2, capsize=0)

errorbar(2,...
    mean(itiByBlk.z(:,3)),...
    std(itiByBlk.z(:,3)),...
    'ko', markerfacecolor='k', linewidth=2, capsize=0)

errorbar(3,...
    mean(itiByBlk.z(:,5)),...
    std(itiByBlk.z(:,5)),...
    'ro', markerfacecolor='r', linewidth=2, capsize=0)
xlim([0.5 3.5])

xticks(1:3)
xticklabels({'Low', 'Mixed', 'High'})

xlabel('Block')
ylabel('Trial initiation time (z-score)')

title({length(ratList), pITIbyBlkPop})

%--------------------------------------------------------------------------
% Regression coefficients
%--------------------------------------------------------------------------
plotThese = 7;
xs = 1:plotThese;

subplot(ht, wd, 4); hold on
plot(xs, betasAll(:, xs)',...
    color=[0.7 0.7 0.7 0.15], linewidth=0.5)
plot(xs, mean(betasAll(:, xs)), 'k', linewidth=2)

plot(xs(pRegressionPop(1:plotThese) < 0.05), 0.01, 'ko')
ylim([-0.1 0.01])

xticks(1:plotThese)
xlabel('Trials back')
ylabel('Regression coefficient')

set(gcf, Position=[496 687 1185 258], renderer='painters')
end