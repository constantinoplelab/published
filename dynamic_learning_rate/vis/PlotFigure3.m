function [] = PlotFigure3(datadir, codedir)
%PlotFigure3 - Plots Figure 3. Must be run after ProcessData_Figure3.
% INPUTS:
%   datadir - Directory where output from ProcessData_Figure3 was saved
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
% To see pipeline for Figure 3, see Process/ProcessData_Figure3.m or 
% Process/ProcessPaperData.m
load([datadir, 'Figure3_Data'], 'BeliefExample',...
    'G_DB1', 'G_PH1', 'G_MS1', 'G_db',...
    'varRatioRat', 'varRatioPH', 'varRatioMS', 'varRatioDB',...
    'dITIbyG_small', 'dITIbyG_large', 'dITIbyG_small_mean_DB',...
    'dITIbyG_small_sem_DB', 'dITIbyG_large_mean_DB',...
    'dITIbyG_large_sem_DB', 'dITIbyG_small_mean_PH',...
    'dITIbyG_small_sem_PH', 'dITIbyG_large_mean_PH', ...
    'dITIbyG_large_sem_PH', 'dITIbyG_small_mean_MS',...
    'dITIbyG_large_mean_MS', 'dITIbyG_small_sem_MS',...
    'dITIbyG_large_sem_MS', 'rpes', 'iexample', 'negBnds', 'posBnds',...
    'G_msByChangePoint', 'BOCD_nTrials',...
    'Bayes_nTrials_binned', 'RunLengthPosterior', 'G_dbByChangePoint',...
    'G_phByChangePoint', 'r', 'p');

%% Plot
nrows = 4;
ncols = 4;

figure; set(gcf, Position=[207 108 1474 837])

%--------------------------------------------------------------------------
% Value examples
%--------------------------------------------------------------------------

subplot(nrows, ncols, 2); hold on
fill([0 40 40 0], [0 0 1 1], 'r', facealpha=0.15, linestyle='none')
fill([40 80 80 40], [0 0 1 1], 'k', facealpha=0.15, linestyle='none')
fill([80 120 120 80], [0 0 1 1], 'b', facealpha=0.15, linestyle='none')

plot(BeliefExample(1,:), 'k')
plot(BeliefExample(2,:), 'r')
plot(BeliefExample(3,:), 'b')

xlim([30 90])

%--------------------------------------------------------------------------
% Value examples
%--------------------------------------------------------------------------
subplot(nrows, ncols, 3); hold on
fill([0 40 40 0], [0 0 1 1], 'r', facealpha=0.15, linestyle='none')
fill([40 80 80 40], [0 0 1 1], 'k', facealpha=0.15, linestyle='none')
fill([80 120 120 80], [0 0 1 1], 'b', facealpha=0.15, linestyle='none')

N = max(G_DB1);

plot(G_DB1/N, 'k')
plot(G_PH1/N, 'r')
plot(G_MS1/N, 'b')

legend('', '', '', 'DB', 'PH', 'MS')
yticks([])
xticks([])

ylabel('Value')
xlabel('Trials')

xlim([30 90])

%--------------------------------------------------------------------------
% Variance
%--------------------------------------------------------------------------
psVar = [ranksum(varRatioRat, varRatioPH),...
    ranksum(varRatioRat, varRatioMS),...
    ranksum(varRatioRat, varRatioDB)];

subplot(nrows, ncols, 4); hold on
errorbar(1:4,...
    mean([varRatioRat, varRatioPH, varRatioMS, varRatioDB]),...
    sem([varRatioRat, varRatioPH, varRatioMS, varRatioDB]),...
    'ko', capsize=0, linewidth=2, markerfacecolor='k')

xticks(1:4)
xticklabels({'Rat', 'PH', 'MS', 'DB'})

xlabel('log-Variance (early/late)')
yline(0, 'k--')

xlim([0.5 4.5])
title(psVar)

%--------------------------------------------------------------------------
% Delta-ITI rat
%--------------------------------------------------------------------------
psdITI = [signrank(dITIbyG_small(:,1), dITIbyG_large(:,1)),...
    signrank(dITIbyG_small(:,2), dITIbyG_large(:,2))];

subplot(nrows, ncols, 5); hold on
shadedErrorBar(1:2,...
    mean(dITIbyG_small, 'omitnan'),...
    sem(dITIbyG_small, 'omitnan'),...
    lineprops={'b-'})

shadedErrorBar(1:2,...
    mean(dITIbyG_large, 'omitnan'),...
    sem(dITIbyG_large, 'omitnan'),...
    lineprops={'r'})

xticks(1:2)
xticklabels({'RPE < 0', 'RPE > 0'})

legend('Low Gain', 'High Gain')
xlim([0.5 2.5])
title('rats')
title({'rats', psdITI})

%--------------------------------------------------------------------------
% model delta-ITI
%--------------------------------------------------------------------------
subplot(nrows, ncols, 6); hold on
shadedErrorBar(1:2,...
    dITIbyG_small_mean_DB,...
    dITIbyG_small_sem_DB, lineprops={'color', 'b'})
shadedErrorBar(1:2,...
    dITIbyG_large_mean_DB,...
    dITIbyG_large_sem_DB, lineprops={'color', 'r'})

xlim([0.5 2.5])
xticks(1:2)
xticklabels({'RPE < 0', 'RPE > 0'})
title('delta-belief')

subplot(nrows, ncols, 7); hold on
shadedErrorBar(1:2,...
    dITIbyG_small_mean_PH,...
    dITIbyG_small_sem_PH, lineprops={'color', 'b'})
shadedErrorBar(1:2,...
    dITIbyG_large_mean_PH,...
    dITIbyG_large_sem_PH, lineprops={'color', 'r'})

xlim([0.5 2.5])
xticks(1:2)
xticklabels({'RPE < 0', 'RPE > 0'})
title('Pearce-Hall')

subplot(nrows, ncols, 8); hold on
shadedErrorBar(1:2, dITIbyG_small_mean_MS, dITIbyG_small_sem_MS,...
    lineprops={'b'})
shadedErrorBar(1:2, dITIbyG_large_mean_MS, dITIbyG_large_sem_MS,...
    lineprops={'r'})

xlim([0.5 2.5])
xticks(1:2)
xticklabels({'RPE < 0', 'RPE > 0'})
title('Mackintosh')

%--------------------------------------------------------------------------
% RPE histogram
%--------------------------------------------------------------------------

subplot(nrows, ncols, 1); hold on
histogram(rpes{iexample}, NumBins = 50);

yl = ylim;

fill([negBnds(65,1) negBnds(65,2) negBnds(65,2) negBnds(65,1)],...
    [yl(1) yl(1) yl(2) yl(2)],...
    'k', facealpha=0.15, linestyle='none')
fill([posBnds(65,1) posBnds(65,2) posBnds(65,2) posBnds(65,1)],...
    [yl(1) yl(1) yl(2) yl(2)],...
    'k', facealpha=0.15, linestyle='none')

xlim([-5 5])

%--------------------------------------------------------------------------
% Change Point posterior
%--------------------------------------------------------------------------

subplot(nrows, ncols, 9); hold on
imagesc(RunLengthPosterior, XData=1:8000)
colorbar
xlim([80-5 80+25])

x1 = 82; plot([x1-0.5 x1+0.5 x1+0.5 x1-0.5 x1-0.5], [0 0 50 50 0], 'w')
x2 = 100; plot([x2-0.5 x2+0.5 x2+0.5 x2-0.5 x2-0.5], [0 0 50 50 0], 'w')
ylim([1 50])

subplot(nrows, ncols, 10); hold on
plot(RunLengthPosterior(:,x1))

xline(x1-80)

xlim([0 40])
yl = ylim;

subplot(nrows, ncols, 11); hold on
plot(RunLengthPosterior(:,x2))

xline(x2-80)

xlim([0 40])
ylim(yl)

%--------------------------------------------------------------------------
% Change Point CDF
%--------------------------------------------------------------------------

subplot(nrows, ncols, 12); hold on
line = cdfplot(BOCD_nTrials);

line.Color='k';

xlim([0 40])

%--------------------------------------------------------------------------
% Comparing BOCD and Bayes
%--------------------------------------------------------------------------
subplot(nrows, ncols, 13); hold on

errorbar(0:10,...
    Bayes_nTrials_binned(1,:),...
    Bayes_nTrials_binned(2,:),...
    'ko', markerfacecolor='k', capsize=0, linewidth=2)

plot([-0.5 10.5], [-0.5, 10.5], 'r--')
xlim([-0.5 10.5])
ylim([-0.5 10.5])

xlabel('BOCD')
ylabel('Bayes')
title([r, p])

%--------------------------------------------------------------------------
% Comparing gain and p(changepoint)
%--------------------------------------------------------------------------

subplot(nrows, ncols, 14); hold on
plot(70:120, RunLengthPosterior(2,70:120))

ylabel('P(Changepoint)')

yyaxis right
plot(70:120, G_db(70:120))

ylim([0 17])
ylabel('Gain')

xlim([70 120])

%--------------------------------------------------------------------------
% Quantifying
%--------------------------------------------------------------------------
subplot(nrows, ncols, 15); hold on
errorbar(1:2, G_dbByChangePoint(1,:), G_dbByChangePoint(2,:),...
    'k.-', capsize=0, linewidth=1, markersize=15)
errorbar(1:2, G_phByChangePoint(1,:), G_phByChangePoint(2,:),...
    'r.-', capsize=0, linewidth=1, markersize=15)
errorbar(1:2, G_msByChangePoint(1,:), G_msByChangePoint(2,:),...
    'b.-', capsize=0, linewidth=1, markersize=15)

yline(1, 'k--')
xlim([0.5 2.5])

xticks(1:2)
xticklabels({'No changepoint', 'Changepoint'})
ylabel('Learning rate gain')

end