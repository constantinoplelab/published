function plotFigure1(savedBehaviorPath)
%Plot figure 1
% YOU MUST RUN processBehaviorData BEFORE RUNNING THIS FUNCTION

%INPUTS:
%   savedBehaviorPath = path to where you have saved the outputs from
%       running processBehaviorData

%load processed behavior data
load([savedBehaviorPath, filesep, 'expert.mat'])
load([savedBehaviorPath, filesep, 'naive.mat'])

%general functions
sem = @(x) std(x,'omitnan')./sqrt(sum(any(~isnan(x), 2)));
setDefaultFigProps

twin = 30; %trial window for wait time dynamics plot

x = 1:5;
xvec = -twin:twin-1;
q = 1:4;
fsize = [10 10 20 8];

%% averages

avgEm = mean(expert.mix, 'omitnan');
semEm= sem(expert.mix);

avgEh = mean(expert.hi, 'omitnan');
semEh = sem(expert.hi);

avgEl = mean(expert.lo, 'omitnan');
semEl = sem(expert.lo);

avgNm = mean(naive.mix, 'omitnan');
semNm = sem(naive.mix);

avgNh = mean(naive.hi, 'omitnan');
semNh = sem(naive.hi);

avgNl = mean(naive.lo, 'omitnan');
semNl = sem(naive.lo);

e_mtol = mean(expert.mtol, 'omitnan');
e_mtol_sem = sem(expert.mtol);

e_mtoh = mean(expert.mtoh, 'omitnan');
e_mtoh_sem = sem(expert.mtoh);

e_ltom = mean(expert.ltom, 'omitnan');
e_ltom_sem = sem(expert.ltom);

e_htom = mean(expert.htom, 'omitnan');
e_htom_sem = sem(expert.htom);

n_mtol = mean(naive.mtol, 'omitnan');
n_mtol_sem = sem(naive.mtol);

n_mtoh = mean(naive.mtoh, 'omitnan');
n_mtoh_sem = sem(naive.mtoh);

n_ltom = mean(naive.ltom, 'omitnan');
n_ltom_sem = sem(naive.ltom);

n_htom = mean(naive.htom, 'omitnan');
n_htom_sem = sem(naive.htom);

e_postLow = mean(expert.postLow, 'omitnan');
e_postLow_sem = sem(expert.postLow);

e_postHigh = mean(expert.postHigh, 'omitnan');
e_postHigh_sem = sem(expert.postHigh);

e_postLowq1 = mean(expert.postLow_q1, 'omitnan');
e_postLowq1_sem = sem(expert.postLow_q1);

e_postHighq1 = mean(expert.postHigh_q1, 'omitnan');
e_postHighq1_sem = sem(expert.postHigh_q1);

n_postLow = mean(naive.postLow, 'omitnan');
n_postLow_sem = sem(naive.postLow);

n_postHigh = mean(naive.postHigh, 'omitnan');
n_postHigh_sem = sem(naive.postHigh);

n_postLowq1 = mean(naive.postLow_q1, 'omitnan');
n_postLowq1_sem = sem(naive.postLow_q1);

n_postHighq1 = mean(naive.postHigh_q1, 'omitnan');
n_postHighq1_sem = sem(naive.postHigh_q1);

%% stats

p_expwt = signrank(expert.lo(:,3), expert.hi(:,3));
p_expL = anova1(expert.postLow(:, 2:end), [], 'off');
p_expH = anova1(expert.postHigh(:, 2:end), [], 'off');
p_expQ1 = (arrayfun(@(x) signrank(expert.postLow_q1(:,x), ...
    expert.postHigh_q1(:,x)), 1:5)).*5;

p_naivewt = signrank(naive.lo(:,3), naive.hi(:,3));
p_naiveL = anova1(naive.postLow(:, 2:end), [], 'off');
p_naiveH = anova1(naive.postHigh(:, 2:end), [], 'off');
p_naiveQ1 = (arrayfun(@(x) signrank(naive.postLow_q1(:,x), ...
    naive.postHigh_q1(:,x)), 1:5)).*5;

%% plot

figure; hold on
tiledlayout(2, 5, 'TileSpacing', 'compact')

%wait time curves
nexttile
shadedErrorBar(x, avgEm, semEm, 'lineprops', {'k', 'linewidth', 1})
hold on
shadedErrorBar(x, avgEh, semEh, 'lineprops', {'-r', 'linewidth', 1})
shadedErrorBar(x, avgEl, semEl, 'lineprops', {'-b', 'linewidth', 1})
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
ylim([10 14.5])
xlabel('Reward offer')
ylabel('Mean wait wime (s)')
ax1 = gca;
ax1.YRuler.TickLabelGapOffset = 1;
% axis square
title('\rm Expert rats')

%wait time dynamics -- expert
nexttile
shadedErrorBar(xvec, e_mtol, e_mtol_sem, 'lineprops', ...
    {'-b', 'linewidth', 1})
hold on
shadedErrorBar(xvec, e_mtoh, e_mtoh_sem, ...
    'lineprops', {'-r', 'linewidth', 1});
ylim([-0.2 0.3])
yl = ylim;
xlim([-25 25]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
yticks([-0.2, 0, 0.2])
xlabel('Trial from block switch');
ylabel('\Delta z-scored wait time');
% axis square
ax2 = gca;
ax2.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off

nexttile
shadedErrorBar(xvec, e_ltom, e_ltom_sem, ...
    'lineprops', {'-b', 'linewidth', 1}); hold on
shadedErrorBar(xvec, e_htom, e_htom_sem, ...
    'lineprops', {'-r', 'linewidth', 1});
ylim([-0.2 0.3])
yl = ylim;
xlim([-25 25]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
yticks([-0.2, 0, 0.2])
xlabel('Trial from block switch');
ylabel('\Delta z-scored wait time');
% axis square
ax3 = gca;
ax3.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off

% quartiles
nexttile
shadedErrorBar(q, e_postLow, e_postLow_sem, 'lineprops', {'color', ...
    [0.1 0.1 0.6], 'linewidth', 1})
hold on
shadedErrorBar(q, e_postHigh, e_postHigh_sem, 'lineprops', {'color', ...
    [0.6 0.1 0.1], 'linewidth', 1})
set(gca, 'TickDir', 'out'); box off;
xlim([0 5]);
ylim([-0.12 0.15])
xticks([1:4])
yticks([-0.1, 0, 0.1])
xlabel('Quartile')
ylabel({'Mean z-scored', 'wait time (s)'})
text(1, 0.14, 'Post-low', 'color', [0.1 0.1 0.6], 'FontSize', 8)
text(1, 0.11, 'Post-high', 'color', [0.6 0.1 0.1], 'FontSize', 8)
ax4 = gca;
ax4.YRuler.TickLabelGapOffset = 1;
% axis square

nexttile
shadedErrorBar(x, e_postLowq1, e_postLowq1_sem, 'lineprops', {'color', ...
    [0.1 0.1 0.6], 'linewidth', 1})
hold on
shadedErrorBar(x, e_postHighq1, e_postHighq1_sem, 'lineprops', {'color', ...
    [0.6 0.1 0.1], 'linewidth', 1})
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
ylim([9.5 15.5])
xlabel('Reward offer')
ylabel('Wait time (s)')
ax5 = gca;
ax5.YRuler.TickLabelGapOffset = 1;
% axis square

nexttile
shadedErrorBar(x, avgNm, semNm, 'lineprops', {'k', 'linewidth', 1})
shadedErrorBar(x, avgNh, semNh, 'lineprops', {'r', 'linewidth', 1})
shadedErrorBar(x, avgNl, semNl, 'lineprops', {'b', 'linewidth', 1})
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
ylim([10 14.5])
xlabel('Reward offer')
ylabel('Mean wait time (s)')
ax6 = gca;
ax6.YRuler.TickLabelGapOffset = 1;
% axis square
title('\rm First 15 sessions')
set(gcf,'renderer','painter')
set(gcf, 'Color', [1 1 1]);

nexttile
shadedErrorBar(xvec, n_mtol, n_mtol_sem, ...
    'lineprops', {'b', 'linewidth', 1});
shadedErrorBar(xvec, n_mtoh, n_mtoh_sem, ...
    'lineprops', {'r', 'linewidth', 1});
ylim([-0.2 0.3])
yl = ylim;
xlim([-25 25]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
yticks([-0.2, 0, 0.2])
xlabel('Trial from block switch');
ylabel('\Delta z-scored wait time');
% axis square
ax7 = gca;
ax7.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off

nexttile
shadedErrorBar(xvec, n_ltom, n_ltom_sem, ...
    'lineprops', {'b', 'linewidth', 1});
shadedErrorBar(xvec, n_htom, n_htom_sem, ...
    'lineprops', {'r', 'linewidth', 1});
ylim([-0.2 0.3])
yl = ylim;
xlim([-25 25]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
yticks([-0.2, 0, 0.2])
xlabel('Trial from block switch');
ylabel('\Delta z-scored wait time');
% axis square
ax8 = gca;
ax8.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off

nexttile
shadedErrorBar(q, n_postLow, n_postLow_sem, 'lineprops', {'color', ...
    [0.1 0.1 0.6], 'linewidth', 1})
hold on
shadedErrorBar(q, n_postHigh, n_postHigh_sem, 'lineprops', {'color', ...
    [0.6 0.1 0.1], 'linewidth', 1})
set(gca, 'TickDir', 'out'); box off;
xlim([0 5]);
ylim([-0.12 0.15])
xticks([1:4])
yticks([-0.1, 0, 0.1])
xlabel('Quartile')
ylabel({'Mean z-scored', 'wait time (s)'})
ax9 = gca;
ax9.YRuler.TickLabelGapOffset = 1;
% axis square

nexttile
shadedErrorBar(x, n_postLowq1, n_postLowq1_sem, 'lineprops', {'color', ...
    [0.1 0.1 0.6], 'linewidth', 1})
hold on
shadedErrorBar(x, n_postHighq1, n_postHighq1_sem, 'lineprops', {'color', ...
    [0.6 0.1 0.1], 'linewidth', 1})
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
ylim([9.5 15.5])
xlabel('Reward offer')
ylabel('Wait time (s)')
ax10 = gca;
ax10.YRuler.TickLabelGapOffset = 1;
% axis square

set(gcf, 'units', 'centimeters', 'position', fsize)
