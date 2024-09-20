function [] = PlotFigure2(datadir, codedir)
%PlotFigure2 - Plots Figure 2. Must be run after ProcessData_Figure2.
% INPUTS:
%   datadir - Directory where output from ProcessData_Figure2 was saved
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
% To see pipeline for Figure 2, see Process/ProcessData_Figure2.m or 
% Process/ProcessPaperData.m
load([datadir, 'Figure2_Data'], 'ratList', 'twin', 'nback',...
    'exampleRatMdl', 'ltom', 'htom', 'betasEarly', 'betasLate',...
    'betasEarlySim', 'betasLateSim',...
    'tauEarly', 'tauLate', 'nLLEarlyEarly', 'nLLEarlyLate',...
    'nLLLateEarly', 'nLLLateLate', 'expParamsEarly', 'expParamsLate',...
    'ratBlk', 'mdlBlk', 'alphaEarly', 'alphaLate', 'Nrats', 'V')

%%
earlyColor =  '#A24289';
lateColor = '#40A158';

f = @(x, params) params(1)*exp(-x/params(2));

exampleRatExpI = 186;

figure
subplot(1, 2, 1); hold on
plot(1:7, betasEarly(exampleRatExpI,1:7),...
    color=earlyColor)
plot(1:7, f(1:7, expParamsEarly(exampleRatExpI,:)),...
    '--', color=earlyColor)
title({ratList{exampleRatExpI}, expParamsEarly(exampleRatExpI,:)})

ylim([-0.1 0.005])

subplot(1, 2, 2); hold on
plot(1:7, betasLate(exampleRatExpI,1:7),...
    color=lateColor)
plot(1:7, f(1:7, expParamsLate(exampleRatExpI,:)),...
    '--', color=lateColor)
title({ratList{exampleRatExpI}, expParamsLate(exampleRatExpI,:)})
ylim([-0.1 0.005])

set(gcf, Position=[1236 771 635 198])

%%  
clear l
ht = 3;
wd = 3;

earlyColor =  '#A24289';
lateColor = '#40A158';

earlyColorRGB =...
    [0.6352941176470588, 0.25882352941176473, 0.5372549019607843];
lateColorRGB =...
    [0.25098039215686274, 0.6313725490196078, 0.34509803921568627];

figure

%--------------------------------------------------------------------------
% Block dynamics
%--------------------------------------------------------------------------
n = 12;
subplot(ht, wd, 1); hold on
fill([0 n n 0], [-0.2 -0.2 0.2 0.2],...
    earlyColorRGB, linestyle='none', facealpha=0.1)
fill([n 40 40 n], [-0.2 -0.2 0.2 0.2],...
    'k', linestyle='none', facealpha=0.1)

l(1) = shadedErrorBar(-twin:twin,...
    mean(ltom), std(ltom)./Nrats, lineprops={'-b'});
l(2) = shadedErrorBar(-twin:twin,...
    mean(htom), std(htom)./Nrats, lineprops={'-r'});
xlim([-20, 40])

%--------------------------------------------------------------------------
% regerssion - early
%--------------------------------------------------------------------------

xs = 1:7;

psEarly = arrayfun(@(i) signrank(betasEarly(:,i)), xs);

subplot(ht, wd, 2); hold on
plot(xs, betasEarly(:,xs), color=[0.7 0.7 0.7 0.1], linewidth=0.25)
plot(xs,...
    mean(betasEarly(:, xs), 'omitnan'),...
    color=earlyColor, linewidth=1);

plot(xs(psEarly<0.05), 0.0025, 'o', color=earlyColor)

xticks(xs)

xlim([xs(1)-0.5, xs(end)])
ylim([-0.2000, 0.0025])
title({'Early', length(ratList)})

%--------------------------------------------------------------------------
% regression - late
%--------------------------------------------------------------------------

psLate = arrayfun(@(i) signrank(betasLate(:,i)), xs);

subplot(ht, wd, 3); hold on
plot(xs, betasLate(:, xs), color=[0.7 0.7 0.7 0.15], linewidth=0.25)
plot(xs,...
    mean(betasLate(:, xs), 'omitnan'),...
    color=lateColor, linewidth=1);

plot(xs(psLate<0.05), 0.001, 'o', color=lateColor)

ylim([-0.2000, 0.0025])
xticks(xs)
xlim([xs(1)-0.5, xs(end)])
title({'Late', length(ratList)})

%--------------------------------------------------------------------------
% taus - early vs. late
%--------------------------------------------------------------------------

pTau = signrank(tauEarly, tauLate, tail='left');
binEdges = (-7:0.5:7) + 0.25;

subplot(ht, wd, 4); hold on
histogram(tauEarly-tauLate, binEdges=binEdges)

xline(0, 'k--', alpha=1, linewidth=1)
xline(mean(tauEarly-tauLate, 'omitnan'), 'r-', alpha=1, linewidth=1)

xlabel('tauEarly')
xlabel('tauLate')

title(pTau)


%--------------------------------------------------------------------------
% model 
%--------------------------------------------------------------------------
subplot(ht, wd, 5); hold on

fill([0 15 15 0], [0 0 1 1], 'r', facealpha = 0.1, linestyle='none')
fill([15 55 55 15], [0 0 1 1], 'k', facealpha = 0.1, linestyle='none')
fill([55 70 70 55], [0 0 1 1], 'b', facealpha = 0.1, linestyle='none')

plot(V, 'k')
plot([15 25], [1 1], color=earlyColor)
plot([45 55], [1 1], color=lateColor)

xlim([0 70])

xticks([])
yticks([])

xlabel('Trials')
ylabel('Value of Environment')
%--------------------------------------------------------------------------
% model vs. data
%--------------------------------------------------------------------------
subplot(ht, wd, 6); hold on
errorbar(0.85, ratBlk(1, 3), ratBlk(2, 3),...
    'bo', markerfacecolor='b', linewidth=2, capsize=5)
errorbar(1.15, mdlBlk(1, 3), ratBlk(2, 3),...
    'bo', markerfacecolor='b', linewidth=2, capsize=0)

errorbar(1.85, ratBlk(1, 1), ratBlk(2, 1),...
    'ko', markerfacecolor='k', linewidth=2, capsize=5)
errorbar(2.15, mdlBlk(1, 1), ratBlk(2, 1),...
    'ko', markerfacecolor='k', linewidth=2, capsize=0)

errorbar(2.85, ratBlk(1, 2), ratBlk(2, 2),...
    'ro', markerfacecolor='r', linewidth=2, capsize=5)
errorbar(3.15, mdlBlk(1, 2), ratBlk(2, 2),...
    'ro', markerfacecolor='r', linewidth=2, capsize=0)

title([exampleRatMdl ' Test Data (Cap)'])
xlim([0.5 3.5])
ylim([0.4 1])
xticks(1:3)
xticklabels({'low', 'mixed', 'high'})
xlabel('Block')

%--------------------------------------------------------------------------
% alpha
%--------------------------------------------------------------------------

subplot(ht, wd, 7); hold on

plot([1 2], [alphaEarly, alphaLate],...
    color=[0.7 0.7 0.7 0.15], linewidth=0.5)
plot([1 2],...
    mean([alphaEarly, alphaLate]),...
    'k', linewidth=2)
xlim([0.75 2.25])
ylim([0 1])

xticks(1:2)
xticklabels({'early', 'late'})
ylabel('learning rate')
title(signrank(alphaEarly, alphaLate))

%--------------------------------------------------------------------------
% simulated regression
%--------------------------------------------------------------------------

psEarlySim = arrayfun(@(i) signrank(betasEarlySim(:,i)), 1:nback);
psLateSim = arrayfun(@(i) signrank(betasLateSim(:,i)), 1:nback);

subplot(ht, wd, 8); hold on

shadedErrorBar(xs,...
    mean(betasEarlySim(:,xs)),...
    sem(betasEarlySim(:,xs)),...
    lineprops={'color', earlyColor}); 
shadedErrorBar(xs,...
    mean(betasLateSim(:,xs)),...
    sem(betasLateSim(:,xs)),...
    lineprops={'color', lateColor}); 

plot(find(psEarlySim<0.05), 0.005, 'o', color=earlyColor)
plot(find(psLateSim<0.05), 0.0025, 'o', color=lateColor)

xticks(xs)
xlim([xs(1)-0.5 xs(end)+0.5])
ylim([-0.11, 0.005])
title('regression for model fits')
%--------------------------------------------------------------------------
% mdl comparsion
%--------------------------------------------------------------------------
earlyTest = nLLEarlyEarly-nLLLateEarly;
lateTest = nLLEarlyLate-nLLLateLate;

pEarly = signrank(earlyTest);
pLate = signrank(lateTest);

subplot(ht, wd, 9); hold on

violinplot([earlyTest lateTest], {'Early', 'Late'},...
    'ShowBox', false, 'ShowWhiskers', false,...
    'ShowMedian', true, 'MarkerSize', 5,...
    'ViolinColor', [0.6429 0.2619 0.5437; 0.2540 0.6389 0.3492]);

yline(0, 'k--', linewidth=1, alpha=1)

xticks(1:2)

xticklabels({'Early Train', 'Late Train'})
ylabel('Early - Late Test nLL')
% ylim([-1 1])

title([earlyTest lateTest])


arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), l)
set(gcf, Position=[559 81 954 785], renderer='painters')
end