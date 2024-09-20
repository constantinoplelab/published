function [] = PlotFigure4(datadir, codedir)
%PlotFigure4 - Plots Figure 4. Must be run after ProcessData_Figure4.
% INPUTS:
%   datadir - Directory where output from ProcessData_Figure4 was saved
%   codedir - Directory of code (e.g., published/dynamic_learning_rate/)

%% Set Path and Load Data

s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codedir, 'IgnoreCase', ispc);
if ~onPath % only add code dir to path if it already isn't
    addpath(genpath(codedir))
end

% Load pre-processed data
% To see pipeline for Figure 4, see Process/ProcessData_Figure4.m or 
% Process/ProcessPaperData.m
load([datadir 'Figure4_Data'], 'CPInByRewMean', 'CPInByRewSEM',...
    'betas_RewardHistory_All', ...
    'betas_RewardHistory_Early', 'betas_RewardHistory_Late', ...
    'CPInByRewMean_All', 'CPInByRewSEM_All', 'CPInByRewAUC', ...
    'CPIn20ByBlkMean_All', 'CPIn20ByBlkSEM_All', 'CPInBy20BlkAUC', ...
    'xsAUC', 'AUCbyRPE_LateMean', 'betas_RPE_EarlyPos', ...
    'betas_RPE_EarlyNeg', 'betas_RPE_LatePos', 'betas_RPE_LateNeg',...
    'RPEbySess_early', 'AUCbySess_early', 'RPEbySess_late', 'AUCbySess_late',...
    'betas_RPE_EarlyPos_All', 'betas_RPE_LatePos_All',...
    'betas_RPE_EarlyNeg_All', 'betas_RPE_LateNeg_All', 'xsDelay', ...
    'AUCBaselineCorrected','sideOnByDelayMean', 'sideOnByDelaySEM',...
    'sideOffByDelayMean', 'sideOffByDelaySEM', 'ratList', 'timesCPIn',...
    'timesSideOn', 'timesSideOff', 'nback', 'bins')

%%
clear l

earlyColor =  '#A24289'; 
lateColor = '#40A158';

ht = 4;
wd = 3;
cmatBlk = [0 0 0; 1 0 0; 0 0 1];
cmat = generate_color_gradient([0 0 1], [1 0 0], 5);

exampleRat = 'G024';
iRat = strcmp(ratList, exampleRat);

figure
set(gcf, Position=[679 42 1002 954], renderer='painter')


subplot(ht, wd, 1); hold on
xs = 1:6;
is = 3:8;

l(1) = shadedErrorBar(xs,...
    mean(betas_RewardHistory_All(:,is)),...
    sem(betas_RewardHistory_All(:,is)),...
    lineprops={'color', 'k'});
l(end+1) = shadedErrorBar(xs,...
    mean(betas_RewardHistory_Early(:,is)),...
    sem(betas_RewardHistory_Early(:,is)),...
    lineprops={'color', earlyColor});
l(end+1) = shadedErrorBar(xs,...
    mean(betas_RewardHistory_Late(:,is)),...
    sem(betas_RewardHistory_Late(:,is)),...
    lineprops={'color', lateColor});
xticks(xs)
xlim([0.5 6.5])

%--------------------------------------------------------------------------
% Plot average CPIn response for one rat
%--------------------------------------------------------------------------
subplot(ht, 5, 3); hold on
fill([0 0.5 0.5 0], [-1.3 -1.3 3 3],...
    'k', linestyle='none', facealpha=0.1)
for r = 1:5
    l(end+1) =...
        shadedErrorBar(timesCPIn,...
        CPInByRewMean{iRat}(r,:),...
        CPInByRewSEM{iRat}(r,:),...
        lineprops={'color', cmat(r,:), 'linewidth', 1});
end

xlim([-0.5 1.15])
ylim([-1.3 3])

xlabel('Time from offer cue (s)')
ylabel('DA (z-score)')
title(exampleRat)

%--------------------------------------------------------------------------
% Plot average CPIn response for all rats
%--------------------------------------------------------------------------
subplot(ht, 5, 4); hold on
fill([0 0.5 0.5 0], [-1 -1 2.1 2.1],...
    'k', facealpha=0.1, linestyle='none')

for r = 1:5
    l(end+1) = shadedErrorBar(timesCPIn,...
        CPInByRewMean_All(r,:), CPInByRewSEM_All(r,:),...
        lineprops={'color', cmat(r,:), 'linewidth', 1});
end

xlim([-0.5 1.15])
ylim([-0.7 2])

xlabel('Time from offer cue (s)')
title(sprintf('N = %d', length(ratList)))

%--------------------------------------------------------------------------
% CPIn response for 20 uL by block
%--------------------------------------------------------------------------
subplot(ht, 5, 5); hold on
fill([0 0.5 0.5 0], [-0.7 -0.7 0.6 0.6],...
    'k', facealpha=0.1, linestyle='none')

for blk = 2:3
    l(end+1) = shadedErrorBar(timesCPIn,...
        CPIn20ByBlkMean_All(blk,:),...
        CPIn20ByBlkSEM_All(blk,:),...
        lineprops={'color', cmatBlk(blk,:), 'linewidth', 1});

end
xlim([-0.5 1.15])
ylim([-0.7 0.6])

xlabel('Time from offer cue (s)')
ylabel('DA (z-score)')
title(sprintf('N = %d', length(ratList)))

%--------------------------------------------------------------------------
% CPIn AUC by reward in mixed
%--------------------------------------------------------------------------

subplot(ht, 6, 7); hold on
plot(1:5, CPInByRewAUC, color = [0.7 0.7 0.7 0.25])

for r = 1:5
    errorbar(r, mean(CPInByRewAUC(:, r)), sem(CPInByRewAUC(:, r)),... 
        '.', color=cmat(r,:), markersize=15, capsize=0)
end

yline(0, 'k--', linewidth=1, alpha=1)
xticks(1:5)
xticklabels({'5', '10', '20', '40', '80'})
xlim([0.5 5.5])

%--------------------------------------------------------------------------
% CPIn AUC for 20 uL by block
%--------------------------------------------------------------------------
subplot(ht, 6, 8); hold on
plot(2:3, CPInBy20BlkAUC(:, 3:-1:2), color = [0.7 0.7 0.7 0.25])

errorbar(2,...
    mean(CPInBy20BlkAUC(:,3)),...
    sem(CPInBy20BlkAUC(:, 3)),...
    '.', color='b', markersize=15, capsize=0)
errorbar(3,...
    mean(CPInBy20BlkAUC(:,2)),...
    sem(CPInBy20BlkAUC(:, 2)),...
    '.', color='r', markersize=15, capsize=0)

xticks(2:3)
xticklabels({'Low', 'High'})
yline(0, 'k--', linewidth=1, alpha=1)

title(signrank(CPInBy20BlkAUC(:,2), CPInBy20BlkAUC(:,3)))

xlim([1.5, 3.5])

%--------------------------------------------------------------------------
% Regression
%--------------------------------------------------------------------------
xs = 0:6;
usethese = 2:8;
regressionAllPs = arrayfun(@(i) signrank(betas_RewardHistory_All(:,i)), 2:nback+1);

subplot(ht, wd, 5); hold on
l(end+1) = shadedErrorBar(xs,...
    mean(betas_RewardHistory_All(:, usethese)),...
    sem(betas_RewardHistory_All(:, usethese)),...
    lineprops={'color', 'k', 'linewidth', 1});

plot(find(regressionAllPs<0.05)-1, 0, 'ko')

xticks(xs)
xlim([-0.5 xs(end)+0.5])
ylim([-0.13, 0.35])

xlabel('Trials Back')
ylabel('Regression coefficient')
title(sprintf('N = %d', length(ratList)))

%--------------------------------------------------------------------------
% AUC by model RPE
%--------------------------------------------------------------------------

subplot(ht, wd, 6)

l(end+1) = shadedErrorBar(xsAUC,...
    mean(AUCbyRPE_LateMean, 'omitnan'),...
    sem(AUCbyRPE_LateMean, 'omitnan'),...
    lineprops={'k', 'linewidth', 1});

yline(0, 'k--', alpha=1)

% ylim([-0.19, 0.8])

xlabel('DA AUC')
ylabel('Late Model RPE')

%--------------------------------------------------------------------------
% dopamine by belief - regression
%--------------------------------------------------------------------------
earlyLateBetaPs =...
    arrayfun(@(i) signrank(betas_RewardHistory_Early(:,i),...
    betas_RewardHistory_Late(:,i),...
    tail = 'right'), 2:8);

subplot(ht, wd, 7); hold on
l(end+1) = shadedErrorBar(0:6,...
    mean(betas_RewardHistory_Early(:, 2:8), 'omitnan'),...
    sem(betas_RewardHistory_Early(:, 2:8), 'omitnan'),...
    lineprops={'color', earlyColor, 'linewidth', 1});
l(end+1) = shadedErrorBar(0:6,...
    mean(betas_RewardHistory_Late(:, 2:8), 'omitnan'),...
    sem(betas_RewardHistory_Late(:, 2:8), 'omitnan'),...
    lineprops={'color', lateColor, 'linewidth', 1});

text(1, 0,...
    arrayfun(@num2str, earlyLateBetaPs, 'UniformOutput', false))

xlim([-0.5 6.5])

xlabel('Trials Back')
ylabel('Regression coefficient')
title(sprintf('N = %d', length(ratList)))

%--------------------------------------------------------------------------
% dopamine by belief - AUC
%--------------------------------------------------------------------------
subplot(4, 6, 15); hold on
plot(RPEbySess_early{3}{11}, AUCbySess_early{3}{11}, 'k.')

xs = xlim;

plot([xs(1) 0],...
    betas_RPE_EarlyNeg{3}(11, 1) +...
    betas_RPE_EarlyNeg{3}(11, 2)*[xs(1) 0],...
    '--', color=earlyColor)
plot([0 xs(end)],...
    betas_RPE_EarlyPos{3}(11, 1) +...
    betas_RPE_EarlyPos{3}(11, 2)*[0 xs(end)],...
    '--', color=earlyColor)

title([betas_RPE_EarlyNeg{3}(11, 2), betas_RPE_EarlyPos{3}(11, 2)])

subplot(4, 6, 16); hold on
plot(RPEbySess_late{3}{11}, AUCbySess_late{3}{11}, 'k.')

xs = xlim;

plot([xs(1) 0],...
    betas_RPE_LateNeg{3}(11, 1) +...
    betas_RPE_LateNeg{3}(11, 2)*[xs(1) 0],...
    '--', color=lateColor)
plot([0 xs(end)],...
    betas_RPE_LatePos{3}(11, 1) +...
    betas_RPE_LatePos{3}(11, 2)*[0 xs(end)],...
    '--', color=lateColor)

title([betas_RPE_LateNeg{3}(11, 2), betas_RPE_LatePos{3}(11, 2)])

%--------------------------------------------------------------------------
% dopamine by belief - AUC
%--------------------------------------------------------------------------
binedges = linspace(-1.5, 1.5, 30);

subplot(ht, wd, 9); hold on
h1 = histogram(betas_RPE_EarlyNeg_All(:,2) - betas_RPE_LateNeg_All(:,2),...
    facecolor='b', binedges=binedges);
histogram(betas_RPE_EarlyPos_All(:,2) - betas_RPE_LatePos_All(:,2),...
    facecolor='r', binedges=binedges)


xline(0, 'k--', linewidth=1, alpha=1)
xline(mean(betas_RPE_EarlyPos_All(:,2) -...
    betas_RPE_LatePos_All(:,2), 'omitnan'), 'r-')
xline(mean(betas_RPE_EarlyNeg_All(:,2) -...
    betas_RPE_LateNeg_All(:,2), 'omitnan'), 'b-')

xlim([-1.5 1.5])
title({signrank(betas_RPE_EarlyPos_All(:,2),...
    betas_RPE_LatePos_All(:,2), tail='right')
signrank(betas_RPE_EarlyNeg_All(:,2),...
betas_RPE_LateNeg_All(:,2), tail='right')})

xlabel('Slope (early-late)')
ylabel('N (sessions)')

%--------------------------------------------------------------------------
% side on by delay
%--------------------------------------------------------------------------
cmatDelay = generate_color_gradient([0.6, 0.6, 0.6], [0 0 0], 6);

subplot(ht, wd, 10); hold on
for r = 1:6
    l(end+1) = shadedErrorBar(timesSideOn,...
        sideOnByDelayMean(r,:),...
        sideOnByDelaySEM(r,:),...
        lineprops = {'color', cmatDelay(r,:), 'linewidth', 1});

end
xline(0, 'k--', linewidth=1, alpha=1)
xlim([-0.5, 10])
ylim([-0.09, 0.925])
legend(arrayfun(@(f)...
    sprintf('%4.2f < delay < %4.2f', bins(f), bins(f+1)),...
    1:6, UniformOutput=false))

title('SideOn')

%--------------------------------------------------------------------------
% side off by delay
%--------------------------------------------------------------------------

subplot(ht, wd, 11); hold on
for r = 2:6
    l(end+1) = shadedErrorBar(timesSideOff,...
        sideOffByDelayMean(r,:),...
        sideOffByDelaySEM(r,:),...
        lineprops = {'color', cmatDelay(r,:), 'linewidth', 1});

end
xlim([-1 1])
ylim([-0.25 1.7])
xline(0, 'k--', linewidth=1, alpha=1)


subplot(ht, wd, 12); hold on
for dd = 1:6
    errorbar(xsDelay(dd),...
        mean(AUCBaselineCorrected(:,dd)),...
        sem(AUCBaselineCorrected(:,dd)),...
        '.', color=cmatDelay(dd,:), markersize=15, capsize=0, linewidth=1)
end
plot(xsDelay(1:6), AUCBaselineCorrected(:, 1:6), color=[0.7 0.7 0.7 0.15])

xlabel('Binned Delay (s)')
ylabel('Baseline-corrected AUC')

xlim([0.25, 8.75])
ylim([-0.2 1.5])


arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), l)
end
