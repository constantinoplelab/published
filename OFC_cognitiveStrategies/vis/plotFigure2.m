function plotFigure1(behaviorDataPath)


twin = 30; %trial window for wait time dynamics plot

x = 1:5;
xvec = -twin:twin-1;
q = 1:4;
fsize = [10 10 20 10];

%% Load processed behavior data
load([behaviorDataPath, filesep, 'infModel.mat'])
load([behaviorDataPath, filesep, 'divModel.mat'])


inf_mtol = mean(infModel.mtol, 'omitnan');
inf_mtol_sem = sem(infModel.mtol);

inf_mtoh = mean(infModel.mtoh, 'omitnan');
inf_mtoh_sem = sem(infModel.mtoh);

inf_ltom = mean(infModel.ltom, 'omitnan');
inf_ltom_sem = sem(infModel.ltom);

inf_htom = mean(infModel.htom, 'omitnan');
inf_htom_sem = sem(infModel.htom);

div_mtol = mean(divModel.mtol, 'omitnan');
div_mtol_sem = sem(divModel.mtol);

div_mtoh = mean(divModel.mtoh, 'omitnan');
div_mtoh_sem = sem(divModel.mtoh);

div_ltom = mean(divModel.ltom, 'omitnan');
div_ltom_sem = sem(divModel.ltom);

div_htom = mean(divModel.htom, 'omitnan');
div_htom_sem = sem(divModel.htom);

inf_postLow = mean(infModel.postLow, 'omitnan');
inf_postLow_sem = sem(infModel.postLow);

inf_postHigh = mean(infModel.postHigh, 'omitnan');
inf_postHigh_sem = sem(infModel.postHigh);

inf_postLowq1 = mean(infModel.postLow_q1, 'omitnan');
inf_postLowq1_sem = sem(infModel.postLow_q1);

inf_postHighq1 = mean(infModel.postHigh_q1, 'omitnan');
inf_postHighq1_sem = sem(infModel.postHigh_q1);

div_postLow = mean(divModel.postLow, 'omitnan');
div_postLow_sem = sem(divModel.postLow);

div_postHigh = mean(divModel.postHigh, 'omitnan');
div_postHigh_sem = sem(divModel.postHigh);

div_postLowq1 = mean(divModel.postLow_q1, 'omitnan');
div_postLowq1_sem = sem(divModel.postLow_q1);

div_postHighq1 = mean(divModel.postHigh_q1, 'omitnan');
div_postHighq1_sem = sem(divModel.postHigh_q1);

%% plot

figure; hold on
tiledlayout(2, 4, 'TileSpacing', 'compact')

%wait time dynamics -- expert
nexttile
shadedErrorBar(xvec, inf_mtol, e_mtol_sem, 'lineprops', ...
    {'-b', 'linewidth', 1})
hold on
shadedErrorBar(xvec, inf_mtoh, inf_mtoh_sem, ...
    'lineprops', {'-r', 'linewidth', 1});
ylim([-0.2 0.3])
yl = ylim;
xlim([-25 25]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('\Delta z-scored wait time');
% axis square
ax1 = gca;
ax1.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off

nexttile
shadedErrorBar(xvec, inf_ltom, inf_ltom_sem, ...
    'lineprops', {'-b', 'linewidth', 1}); hold on
shadedErrorBar(xvec, inf_htom, inf_htom_sem, ...
    'lineprops', {'-r', 'linewidth', 1});
ylim([-0.2 0.3])
yl = ylim;
xlim([-25 25]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('\Delta z-scored wait time');
% axis square
ax2 = gca;
ax2.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off

% quartiles
nexttile
shadedErrorBar(q, inf_postLow, inf_postLow_sem, 'lineprops', {'color', ...
    [0.1 0.1 0.6], 'linewidth', 1})
hold on
shadedErrorBar(q, inf_postHigh, inf_postHigh_sem, 'lineprops', {'color', ...
    [0.6 0.1 0.1], 'linewidth', 1})
set(gca, 'TickDir', 'out'); box off;
xlim([0 5]);
ylim([-0.12 0.15])
xlabel('Quartile')
ylabel({'Mean z-scored', 'wait wime (s)'})
text(1, 0.14, 'Post-low', 'color', [0.1 0.1 0.6], 'FontSize', 8)
text(1, 0.11, 'Post-high', 'color', [0.6 0.1 0.1], 'FontSize', 8)
ax3 = gca;
ax3.YRuler.TickLabelGapOffset = 1;
% axis square

nexttile
shadedErrorBar(x, inf_postLowq1, inf_postLowq1_sem, 'lineprops', {'color', ...
    [0.1 0.1 0.6], 'linewidth', 1})
hold on
shadedErrorBar(x, inf_postHighq1, inf_postHighq1_sem, 'lineprops', {'color', ...
    [0.6 0.1 0.1], 'linewidth', 1})
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
% ylim([10 14.5])
xlabel('Reward offer')
ylabel('Wait time (s)')
ax4 = gca;
ax4.YRuler.TickLabelGapOffset = 1;
% axis square

nexttile
shadedErrorBar(xvec, div_mtol, div_mtol_sem, ...
    'lineprops', {'b', 'linewidth', 1});
shadedErrorBar(xvec, div_mtoh, div_mtoh_sem, ...
    'lineprops', {'r', 'linewidth', 1});
ylim([-0.2 0.3])
yl = ylim;
xlim([-25 25]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('\Delta z-scored wait time');
% axis square
ax5 = gca;
ax5.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off

nexttile
shadedErrorBar(xvec, div_ltom, div_ltom_sem, ...
    'lineprops', {'b', 'linewidth', 1});
shadedErrorBar(xvec, div_htom, div_htom_sem, ...
    'lineprops', {'r', 'linewidth', 1});
ylim([-0.2 0.3])
yl = ylim;
xlim([-25 25]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('\Delta z-scored wait time');
% axis square
ax6 = gca;
ax6.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off

nexttile
shadedErrorBar(q, div_postLow, div_postLow_sem, 'lineprops', {'color', ...
    [0.1 0.1 0.6], 'linewidth', 1})
hold on
shadedErrorBar(q, div_postHigh, div_postHigh_sem, 'lineprops', {'color', ...
    [0.6 0.1 0.1], 'linewidth', 1})
set(gca, 'TickDir', 'out'); box off;
xlim([0 5]);
ylim([-0.12 0.15])
xlabel('Quartile')
ylabel({'Mean z-scored', 'wait wime (s)'})
text(1, 0.14, 'Post-low', 'color', [0.1 0.1 0.6], 'FontSize', 8)
text(1, 0.11, 'Post-high', 'color', [0.6 0.1 0.1], 'FontSize', 8)
ax7 = gca;
ax7.YRuler.TickLabelGapOffset = 1;
% axis square

nexttile
shadedErrorBar(x, div_postLowq1, div_postLowq1_sem, 'lineprops', {'color', ...
    [0.1 0.1 0.6], 'linewidth', 1})
hold on
shadedErrorBar(x, div_postHighq1, div_postHighq1_sem, 'lineprops', {'color', ...
    [0.6 0.1 0.1], 'linewidth', 1})
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
% ylim([10 14.5])
xlabel('Reward offer')
ylabel('Wait time (s)')
ax8 = gca;
ax8.YRuler.TickLabelGapOffset = 1;
% axis square

set(gcf, 'units', 'centimeters', 'position', fsize)
