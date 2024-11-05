function PlotFigure4(datadir, codedir)
%PlotFigure4 - Plots Figure 4. 
% INPUTS:
%   datadir - Local directory where ProcessData_Figure4.mat was saved after running ProcessData_Figure4
%   codedir - Local directory where code (e.g., published/EstrousRPE/) was saved

%% Set up paths
s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codedir, 'IgnoreCase', ispc);
if ~onPath % only add code dir to path if it already isn't
    addpath(genpath(codedir))
end

%% Load data
load([datadir, 'ProcessData_Figure4'], 'NAcc_ratlist', 'photITIs', 'ITIbins',...
    'da_rats', 'AUC_overrats', 'P_byrat', 'Rcorr_byrat', 'bins', 'T',...
    'ITI_legend', 'LatOpts', 'avgITIbyblock_opto', 'errITIbyblock_opto',...
    'avgITIbyblock_nonopto', 'errITIbyblock_nonopto', 'nonopto_norm',...
    'opto_norm', 'isoptorat', 'ratTrial');

%% PLOT %%
figure;
set(gcf, 'Color', [1 1 1], 'Units', 'Inches', 'Position', [0, 0, 17, 9],...
    'PaperUnits', 'Inches', 'PaperSize', [9, 5], 'renderer','painters')

%set generally used variables
initiationcolors = {'#A2142F','#D95319','#EDB120','#FFFF00','#77AC30',...
    '#0072BD','#7E2F8E','k',[0 0 0 0.15]}; %,'#4DBEEE'
optocolors = {'k','#3AAEAD'};
blocknames = {'Low', 'High'};

%--------------------------------------------------------------------------
%4a. DA trace at offer cue by following initiation time bin. No
%baseline-correction, mixed block and violation trials only.
%--------------------------------------------------------------------------
%initiation time histogram inset
nexttile
histogram(log2(photITIs), binwidth=0.05, facecolor='k', edgecolor='k'); hold on
for b=2:length(ITIbins)
    xline(log2(ITIbins(b)), color=initiationcolors{b-1}, linewidth=0.5)
end
yticks(0)
xticks(-12:4:12)
xlim([-10 9])
set(gca, 'TickDir', 'out'); box off
grid off; 
xlabel('Initiation times (s)')

nexttile
da_bins_rats = cell(length(bins), 1);
for j = 2:length(bins)
    da_bin_rats = NaN(length(NAcc_ratlist), length(T));
    for rat = 1:length(NAcc_ratlist)
        da_rat = da_rats{rat};
        da_bin_rats(rat, :) = da_rat{1, j};
    end
    da_bins_rats{j} = da_bin_rats;
end
for j = 2:length(bins)
    DA = da_bins_rats{j};
    y = DA;
    err = sem(y);
    thisplot = shadedErrorBar(T, mean(y, 'omitnan'), err,...
        'lineprops', {'Color', initiationcolors{j-1}, 'linewidth', 0.5,...
        'DisplayName', ITI_legend{j}}); hold on
    set(thisplot.edge,'LineWidth',0.001,'LineStyle',':')
end
set(gcf, 'Color', [1 1 1]);
xlim([-0.5 1.25]);
xticks(-1:0.5:1)
set(gca, 'TickDir', 'out'); box off;
xlabel('Time from offer cue (s)');   
ylabel('deltaF/F')
axis square
sgtitle(['N=' num2str(length(NAcc_ratlist))]);
hleg = legend('location', 'best');
title(hleg,' Initiation time (less than or equal to, s)')
line([0 0], ylim, 'Color', [0 0 0], 'LineStyle', '--', HandleVisibility='off');
xline([0.5 0.5], '--k', HandleVisibility='off'); hold on %end of AUC window

%--------------------------------------------------------------------------
%4b. Population average AUC by initiation time bin.
%--------------------------------------------------------------------------
nexttile
scatter(1:length(bins)-1, mean(AUC_overrats(:,2:end), 'omitnan'),...
    100, 'black', 'filled'); hold on
errorbar(1:length(bins)-1, mean(AUC_overrats(:,2:end), 'omitnan'),...
    sem(AUC_overrats(:,2:end)), 'k', capsize=10)
axis square
set(gca, 'TickDir', 'out'); box off
xticks(1:length(bins)-1)
xticklabels(ITI_legend(2:end))
xlim([0.5 length(bins)-0.5])
xlabel('Bin')
alpha(0.3)
ylabel('AUC')
ylim([-0.2 0.27])
yticks(-0.2:0.1:0.3)

%--------------------------------------------------------------------------
%4c. Histogram of pearson correlations by rat. All are significant (p <<< 0.0001).
%--------------------------------------------------------------------------
nexttile
histogram(Rcorr_byrat, facecolor='k', binwidth=0.05)
xlim([-0.4 0.4])
xline(0, '--k')
xline(median(Rcorr_byrat), '--r') %%%ADD MAKE MEDIAN%%%
title(['N=' num2str(length(NAcc_ratlist))])
subtitle(['p=' num2str(P_byrat(2:end)')])
axis square; grid off
set(gca, 'TickDir', 'out'); box off;
set(gcf, 'Color', [1 1 1], 'renderer', 'painters');
xlabel('r')
ylabel('# rats')
yticks(0:1:3)

%--------------------------------------------------------------------------
%4g. Model of initiation time with reinforcement learning model,
% simulation of opto stimulation with additional RPE on 30% of trials
%--------------------------------------------------------------------------
nexttile
for oa=1:2
    LatOpt = LatOpts{oa};
    iti_block_avg = NaN(1,2);
    iti_block_er = NaN(1,2);
    for bl = 2:3
        iti_block_avg(bl) = mean(LatOpt(ratTrial.block==bl), 'omitnan');
        iti_block_er(bl) = std(LatOpt(ratTrial.block==bl), 'omitnan')./...
            sqrt(sum(~isnan(LatOpt(ratTrial.block==bl))));
    end
    y = [iti_block_avg(3) iti_block_avg(2)];
    yer = [iti_block_er(3) iti_block_er(2)];
    x = [1 2];
    plots = NaN(1,2);
    for bl = 1:2
        errorbar(x(bl),y(bl),yer(bl),Color=optocolors{oa},LineWidth=0.5); hold on
        plots(oa, 1) = plot(x(bl),y(bl),'.',Color=optocolors{oa},markersize=30); hold on
    end
end
xlim([.5 2.5]);
set(gca, 'TickDir', 'out'); box off
xticks([1 2])
set(gca, 'xTickLabels', {'Low'; 'High'});
xlabel('Reward block');
ylabel('Initiation time (sec)');
yticks(2:1:4)
ylim([1.5 4.5])
set(gcf, 'Color', [1 1 1]);
set(gca, 'TickDir', 'out'); box off;
axis square

%--------------------------------------------------------------------------
%4h. Example opto rat, G073
%--------------------------------------------------------------------------
nexttile
errorbar(1, avgITIbyblock_nonopto(3), errITIbyblock_nonopto(3),...
    color='k', LineWidth=0.5, capsize=10); hold on
errorbar(2, avgITIbyblock_nonopto(2), errITIbyblock_nonopto(2),...
    color='k', LineWidth=0.5, capsize=10); hold on
plot([1 2], [avgITIbyblock_nonopto(3)...
    avgITIbyblock_nonopto(2)], '.k', markersize=30)
errorbar(1, avgITIbyblock_opto(3), errITIbyblock_opto(3),...
    color=optocolors{2}, LineWidth=0.5, capsize=10); hold on
errorbar(2, avgITIbyblock_opto(2), errITIbyblock_opto(2),...
    color=optocolors{2}, LineWidth=0.5, capsize=10); hold on
plot([1 2], [avgITIbyblock_opto(3) ...
    avgITIbyblock_opto(2)], '.', markersize=30, color=optocolors{2})
xlim([0.5 2.5]);
xticks([1 2])
xticklabels({'Low', 'High'})
xlabel('Reward block')
ylabel('Initiation time (s)')
yticks(0:0.4:4)
grid off
set(gca, 'TickDir', 'out'); box off
axis square 

%--------------------------------------------------------------------------
%4i. Average effect of opto on initiation time for TH-Cre + ChR2 rats
% (opto compared to control sessions)
%--------------------------------------------------------------------------
nexttile
blocks_to_plot = [3 2];

%plot opto rats
controlsess = nonopto_norm(isoptorat,:);
optosess = opto_norm(isoptorat,:);
pvals = NaN(1, 2);
for bl=1:2
    errorbar(bl, mean(controlsess(:, blocks_to_plot(bl)), 'omitnan'),...
        std(controlsess(:, blocks_to_plot(bl)), 'omitnan')./...
        sqrt(sum(~isnan(controlsess(:, blocks_to_plot(bl))))),...
        'Color', 'k', 'LineWidth', 0.5, capsize=10); hold on   
    plot(bl, mean(controlsess(:, blocks_to_plot(bl)), 'omitnan'),...
        '.', MarkerSize=25, Color='k')
    errorbar(bl, mean(optosess(:, blocks_to_plot(bl)), 'omitnan'),...
        std(optosess(:, blocks_to_plot(bl)), 'omitnan')./...
        sqrt(sum(~isnan(optosess(:, blocks_to_plot(bl))))),...
        'Color', '#3AAEAD', 'LineWidth', 0.5, capsize=10); hold on
    plot(bl, mean(optosess(:, blocks_to_plot(bl)), 'omitnan'),...
        '.', MarkerSize=25, Color='#3AAEAD')
    pvals(1, bl) = signrank(controlsess(:, blocks_to_plot(bl)),...
        optosess(:, blocks_to_plot(bl)));
end
low_effectsize = effsize(controlsess(:, 3), optosess(:, 3));
high_effectsize = effsize(controlsess(:, 2), optosess(:, 2));
xlim([0.5 2.5]);
xticks([1 2])
xticklabels(blocknames)
grid off
set(gca, 'TickDir', 'out'); box off
axis square 
title(['TH-ChR2 rats N=' num2str(sum(isoptorat))]) 
subtitle(['p=' num2str(pvals)...
    ', low effect size: ' num2str(low_effectsize)...
    ', high effect size: ' num2str(high_effectsize)])
ylabel('Initiation time (normalized to control sess high)')
yticks(0:1:3)
yline(1, '--k')
ylim([0.7 1.65])

%plot control rats
controlsess = nonopto_norm(~isoptorat,:);
optosess = opto_norm(~isoptorat,:);
plots = NaN(1, 2);
nexttile
pvals = NaN(1, 2);
for bl=1:2
    errorbar(bl, mean(controlsess(:, blocks_to_plot(bl)), 'omitnan'),...
        std(controlsess(:, blocks_to_plot(bl)), 'omitnan')./...
        sqrt(sum(~isnan(controlsess(:, blocks_to_plot(bl))))),...
        'Color', 'k', 'LineWidth', 0.5, capsize=10); hold on   
    plots(1) = plot(bl, mean(controlsess(:, blocks_to_plot(bl)), 'omitnan'),...
        '.', MarkerSize=25, Color='k');
    errorbar(bl, mean(optosess(:, blocks_to_plot(bl)), 'omitnan'),...
        std(optosess(:, blocks_to_plot(bl)), 'omitnan')./...
        sqrt(sum(~isnan(optosess(:, blocks_to_plot(bl))))),...
        'Color', '#3AAEAD', 'LineWidth', 0.5, capsize=10); hold on
    plots(2) = plot(bl, mean(optosess(:, blocks_to_plot(bl)), 'omitnan'),...
        '.', MarkerSize=25, Color='#3AAEAD');
    pvals(1, bl) = signrank(controlsess(:, blocks_to_plot(bl)),...
        optosess(:, blocks_to_plot(bl)));
end
xlim([0.5 2.5]);
xticks([1 2])
xticklabels(blocknames)
grid off
set(gca, 'TickDir', 'out'); box off
axis square 
title(['Sham rats N=' num2str(sum(~isoptorat))])  
ylabel('Initiation time (normalized to control sess high)')
legend(plots, {'control sessions', 'opto sessions'}, 'Location','best')
subtitle(['p=' num2str(pvals)])
yticks(0:1:3)
ylim([0.7 2])
yline(1, '--k')


