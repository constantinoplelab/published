function plotFigure6_naiveBeh(behaviorPath, processedBehaviorPath)
%Plot figure 6 - behavior for the first 15 sessions of training
% YOU MUST RUN processBehaviorData BEFORE RUNNING THIS FUNCTION

%INPUTS:
%   behaviorPath = local path to behavior data downloaded from zenodo
%   savedBehaviorPath = path to where you have saved the outputs from
%       running processBehaviorData

%load processed behavior data
load([processedBehaviorPath, filesep, 'naive.mat'])
load([processedBehaviorPath, filesep, 'divModel.mat'])

%general functions
sem = @(x) std(x,'omitnan')./sqrt(sum(any(~isnan(x), 2)));
setDefaultFigProps

twin = 40; %trial window for wait time dynamics plot
x = 1:5;
xvec = -twin:twin-1;
q = 1:4;
fsize = [10 10 22 10];

%% averages

avgNm = mean(naive.mix, 'omitnan');
semNm = sem(naive.mix);

avgNh = mean(naive.hi, 'omitnan');
semNh = sem(naive.hi);

avgNl = mean(naive.lo, 'omitnan');
semNl = sem(naive.lo);

n_mtol = mean(naive.mtol, 'omitnan');
n_mtol_sem = sem(naive.mtol);

n_mtoh = mean(naive.mtoh, 'omitnan');
n_mtoh_sem = sem(naive.mtoh);

n_ltom = mean(naive.ltom, 'omitnan');
n_ltom_sem = sem(naive.ltom);

n_htom = mean(naive.htom, 'omitnan');
n_htom_sem = sem(naive.htom);

div_mtol = mean(divModel.mtol, 'omitnan');
div_mtol_sem = sem(divModel.mtol);

div_mtoh = mean(divModel.mtoh, 'omitnan');
div_mtoh_sem = sem(divModel.mtoh);

div_ltom = mean(divModel.ltom, 'omitnan');
div_ltom_sem = sem(divModel.ltom);

div_htom = mean(divModel.htom, 'omitnan');
div_htom_sem = sem(divModel.htom);

n_postLow = mean(naive.postLow, 'omitnan');
n_postLow_sem = sem(naive.postLow);

n_postHigh = mean(naive.postHigh, 'omitnan');
n_postHigh_sem = sem(naive.postHigh);

n_postLowq1 = mean(naive.postLow_q1, 'omitnan');
n_postLowq1_sem = sem(naive.postLow_q1);

n_postHighq1 = mean(naive.postHigh_q1, 'omitnan');
n_postHighq1_sem = sem(naive.postHigh_q1);

n_postLowq1_con = mean(naive.postLow_q1Con, 'omitnan');
n_postLowq1_con_sem = sem(naive.postLow_q1Con);

n_postHighq1_con = mean(naive.postHigh_q1Con, 'omitnan');
n_postHighq1_con_sem = sem(naive.postHigh_q1Con);


div_postLow = mean(divModel.postLow, 'omitnan');
div_postLow_sem = sem(divModel.postLow);

div_postHigh = mean(divModel.postHigh, 'omitnan');
div_postHigh_sem = sem(divModel.postHigh);

div_postLowq1 = mean(divModel.postLow_q1, 'omitnan');
div_postLowq1_sem = sem(divModel.postLow_q1);

div_postHighq1 = mean(divModel.postHigh_q1, 'omitnan');
div_postHighq1_sem = sem(divModel.postHigh_q1);

div_postLowq1_con = mean(divModel.postLow_q1Con, 'omitnan');
div_postLowq1_con_sem = sem(divModel.postLow_q1Con);

div_postHighq1_con = mean(divModel.postHigh_q1Con, 'omitnan');
div_postHighq1_con_sem = sem(divModel.postHigh_q1Con);

loAvg = mean(naive.lo(:, 1:3), 2, 'omitnan');
loConAvg = mean(naive.postLow_q1Con(:,1:3), 2, 'omitnan');
hiAvg = mean(naive.hi(:, 3:5), 2, 'omitnan');
hiConAvg = mean(naive.postHigh_q1Con(:,3:5), 2, 'omitnan');

%% stats

p_naivewt = signrank(naive.lo(:,3), naive.hi(:,3));
p_naiveL = anova1(naive.postLow(:, 2:end), [], 'off');
p_naiveH = anova1(naive.postHigh(:, 2:end), [], 'off');
p_naiveQ1 = (arrayfun(@(x) signrank(naive.postLow_q1(:,x), ...
    naive.postHigh_q1(:,x)), 1:5)).*5;
p_naiveQ1Con = signrank(naive.postLow_q1Con(:,3), divModel.postHigh_q1Con(:,3));

p_naiveConLow = signrank(loAvg, loConAvg);
p_naiveConHigh = signrank(hiAvg, hiConAvg);

%% simulate divisively normalized value
load([behaviorPath filesep 'All', filesep, 'ratTrialAll_E001']);

params = [50, 0.15, 60]; %k, alpha, tback for divisive normalization model
rew = convertreward(A.reward);
[normval, ~] = divisiveNorm_fun_2(A, params, 'log');

%% plot

figure; hold on
tiledlayout(2, 6, 'TileSpacing', 'compact')

nexttile
shadedErrorBar(x, avgNm, semNm, 'lineprops', {'k', 'linewidth', 1})
shadedErrorBar(x, avgNh, semNh, 'lineprops', {'r', 'linewidth', 1})
shadedErrorBar(x, avgNl, semNl, 'lineprops', {'b', 'linewidth', 1})
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
ylim([10.8 14.8])
xlabel('Reward offer')
ylabel('Mean wait time (s)')
ax1 = gca;
ax1.YRuler.TickLabelGapOffset = 1;
% axis square
title('\rm First 15 sessions')


nexttile
shadedErrorBar(xvec, n_mtol, n_mtol_sem, ...
    'lineprops', {'b', 'linewidth', 1});
shadedErrorBar(xvec, n_mtoh, n_mtoh_sem, ...
    'lineprops', {'r', 'linewidth', 1});
ylim([-0.2 0.3])
yl = ylim;
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
yticks([-0.2, 0, 0.2])
xlabel('Trial from block switch');
ylabel('Wait time (z-scored)');
% axis square
ax2 = gca;
ax2.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off

nexttile
shadedErrorBar(xvec, n_ltom, n_ltom_sem, ...
    'lineprops', {'b', 'linewidth', 1});
shadedErrorBar(xvec, n_htom, n_htom_sem, ...
    'lineprops', {'r', 'linewidth', 1});
ylim([-0.2 0.3])
yl = ylim;
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
yticks([-0.2, 0, 0.2])
% axis square
ax3 = gca;
ax3.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off

nexttile
shadedErrorBar(q, n_postLow, n_postLow_sem, 'lineprops', {'color', ...
    [0.1 0.1 0.6], 'linewidth', 1})
hold on
shadedErrorBar(q, n_postHigh, n_postHigh_sem, 'lineprops', {'color', ...
    [0.6 0.1 0.1], 'linewidth', 1})
set(gca, 'TickDir', 'out'); box off;
xlim([0 5]);
ylim([-0.1 0.15])
xticks([1:4])
yticks([-0.1, 0, 0.1])
xlabel('Quartile')
ylabel({'Mean z-scored', 'wait time'})
ax4 = gca;
ax4.YRuler.TickLabelGapOffset = 1;
% axis square

nexttile
shadedErrorBar(x, n_postLowq1_con, n_postLowq1_con_sem, 'lineprops', {'color', ...
    [0.1 0.1 0.6], 'linewidth', 1})
hold on
shadedErrorBar(x, n_postHighq1_con, n_postHighq1_con_sem, 'lineprops', {'color', ...
    [0.6 0.1 0.1], 'linewidth', 1})
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
ylim([10.8 14.8])
xlabel('Reward offer')
ylabel('Wait time (s)')
title('Congruent')
ax5 = gca;
ax5.YRuler.TickLabelGapOffset = 1;
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
ylim([10.8 14.8])
xlabel('Reward offer')
ylabel('Wait time (s)')
title('Incong')
ax6 = gca;
ax6.YRuler.TickLabelGapOffset = 1;
% axis square

% Divisively normalized value schematic
t = 493:839; %trials for mixed, low, mixed, high, mixed pattern
colors = [0 0 0; 1 0 0; 0 0 1];

nexttile 
for b = 1:3
    blocks = find(A.block(t)==b); %plot block borders
    plot(t(blocks), 11.8.*ones(1, length(blocks)), '.', 'Color', colors(b,:));
    hold on
end

dotcolors = [0 0 1; 0.07 0.76 0.95; 0.69 0.29 0.76; 0.98 0.45 0.45; 1 0 0];

for r = 1:5
    v = find(rew(t) == r);
    scatter(t(v), normval(t(v)), 10, dotcolors(r, :), 'filled');
    hold on
end
set(gca, 'TickDir', 'out'); box off
xlabel('Trials');
ylabel('Value');
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

nexttile
shadedErrorBar(xvec, div_mtol, div_mtol_sem, ...
    'lineprops', {'b', 'linewidth', 1});
shadedErrorBar(xvec, div_mtoh, div_mtoh_sem, ...
    'lineprops', {'r', 'linewidth', 1});
ylim([-2 2])
yl = ylim;
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('Wait time (z-scored)');
% axis square
ax8 = gca;
ax8.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off

nexttile
shadedErrorBar(xvec, div_ltom, div_ltom_sem, ...
    'lineprops', {'b', 'linewidth', 1});
shadedErrorBar(xvec, div_htom, div_htom_sem, ...
    'lineprops', {'r', 'linewidth', 1});
ylim([-2 2])
yl = ylim;
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
% axis square
ax9 = gca;
ax9.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off

nexttile
shadedErrorBar(q, div_postLow, div_postLow_sem, 'lineprops', {'color', ...
    [0.1 0.1 0.6], 'linewidth', 1})
hold on
shadedErrorBar(q, div_postHigh, div_postHigh_sem, 'lineprops', {'color', ...
    [0.6 0.1 0.1], 'linewidth', 1})
set(gca, 'TickDir', 'out'); box off;
xlim([0 5]);
ylim([-1.5 1.5])
xticks([1:4])
xlabel('Quartile')
ylabel({'Mean z-scored', 'wait time'})
ax10 = gca;
ax10.YRuler.TickLabelGapOffset = 1;
% axis square

nexttile
shadedErrorBar(x, div_postLowq1_con, div_postLowq1_con_sem, 'lineprops', {'color', ...
    [0.1 0.1 0.6], 'linewidth', 1})
hold on
shadedErrorBar(x, div_postHighq1_con, div_postHighq1_con_sem, 'lineprops', {'color', ...
    [0.6 0.1 0.1], 'linewidth', 1})
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
ylim([1 11])
xlabel('Reward offer')
ylabel('Wait time (s)')
title('Congruent')
ax11 = gca;
ax11.YRuler.TickLabelGapOffset = 1;
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
ylim([1 11])
xlabel('Reward offer')
ylabel('Wait time (s)')
title('Incong')
ax12 = gca;
ax12.YRuler.TickLabelGapOffset = 1;
% axis square


set(gcf,'renderer','painter')
set(gcf, 'Color', [1 1 1]);
set(gcf, 'units', 'centimeters', 'position', fsize)
