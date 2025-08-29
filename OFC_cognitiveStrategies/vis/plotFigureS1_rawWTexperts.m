function plotFigureS1_rawWTexperts(processedBehaviorPath)
%Plot figure S1 - raw wait time block transition dynamics for expert rats
% and inference model predictions
% YOU MUST RUN processBehaviorData BEFORE RUNNING THIS FUNCTION

%INPUTS:
%   processedBehaviorPath = path to where you have saved the outputs from
%       running processBehaviorData

%load processed behavior data
load([processedBehaviorPath, filesep, 'expert.mat'])
load([processedBehaviorPath, filesep, 'infModel.mat'])

%general functions
sem = @(x) std(x,'omitnan')./sqrt(sum(any(~isnan(x), 2)));
setDefaultFigProps

twin = 40; %trial window for wait time dynamics plot
xvec = -twin:twin-1;
q = 1:4;
fsize = [10 10 20 8];


%% Averages - congruent volumes exluding 20

e_mtol = mean(expert.mtolRaw, 'omitnan');
e_mtol_sem = sem(expert.mtolRaw);

e_mtoh = mean(expert.mtohRaw, 'omitnan');
e_mtoh_sem = sem(expert.mtohRaw);

e_ltom = mean(expert.ltomRaw, 'omitnan');
e_ltom_sem = sem(expert.ltomRaw);

e_htom = mean(expert.htomRaw, 'omitnan');
e_htom_sem = sem(expert.htomRaw);

inf_mtol = mean(infModel.mtolRaw, 'omitnan');
inf_mtol_sem = sem(infModel.mtolRaw);

inf_mtoh = mean(infModel.mtohRaw, 'omitnan');
inf_mtoh_sem = sem(infModel.mtohRaw);

inf_ltom = mean(infModel.ltomRaw, 'omitnan');
inf_ltom_sem = sem(infModel.ltomRaw);

inf_htom = mean(infModel.htomRaw, 'omitnan');
inf_htom_sem = sem(infModel.htomRaw);


%% Averages - 20 ul only

e_mtol20 = mean(expert.mtol20Raw, 'omitnan');
e_mtol20_sem = sem(expert.mtol20Raw);

e_mtoh20 = mean(expert.mtoh20Raw, 'omitnan');
e_mtoh20_sem = sem(expert.mtoh20Raw);

e_ltom20 = mean(expert.ltom20Raw, 'omitnan');
e_ltom20_sem = sem(expert.ltom20Raw);

e_htom20 = mean(expert.htom20Raw, 'omitnan');
e_htom20_sem = sem(expert.htom20Raw);

inf_mtol20 = mean(infModel.mtol20Raw, 'omitnan');
inf_mtol20_sem = sem(infModel.mtol20Raw);

inf_mtoh20 = mean(infModel.mtoh20Raw, 'omitnan');
inf_mtoh20_sem = sem(infModel.mtoh20Raw);

inf_ltom20 = mean(infModel.ltom20Raw, 'omitnan');
inf_ltom20_sem = sem(infModel.ltom20Raw);

inf_htom20 = mean(infModel.htom20Raw, 'omitnan');
inf_htom20_sem = sem(infModel.htom20Raw);

%% Plot

figure; hold on
tiledlayout(2, 4, 'TileSpacing', 'compact')

%wait time dynamics -- EXPERT, congruent volumes
nexttile
shadedErrorBar(xvec, e_mtol, e_mtol_sem, 'lineprops', ...
    {'-b', 'linewidth', 1})
hold on
shadedErrorBar(xvec, e_mtoh, e_mtoh_sem, ...
    'lineprops', {'-r', 'linewidth', 1});
ylim([8.5 13])
yl = ylim;
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('Wait time (s)');
title('\rm Expert: mixed to adapt')
ax1 = gca;
ax1.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off
axis square

nexttile
shadedErrorBar(xvec, e_ltom, e_ltom_sem, ...
    'lineprops', {'-b', 'linewidth', 1}); hold on
shadedErrorBar(xvec, e_htom, e_htom_sem, ...
    'lineprops', {'-r', 'linewidth', 1});
ylim([8.5 13])
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('Wait time (s)');
title('\rm Expert: adapt to mixed')
ax2 = gca;
ax2.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off
axis square

%wait time dynamics -- EXPERT, 20ul only
nexttile
shadedErrorBar(xvec, e_mtol20, e_mtol20_sem, 'lineprops', ...
    {'-b', 'linewidth', 1})
hold on
shadedErrorBar(xvec, e_mtoh20, e_mtoh20_sem, ...
    'lineprops', {'-r', 'linewidth', 1});
ylim([8.5 13])
yl = ylim;
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('Wait time (s)');
title('\rm Expert 20: mixed to adapt')
ax1 = gca;
ax1.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off
axis square

nexttile
shadedErrorBar(xvec, e_ltom20, e_ltom20_sem, ...
    'lineprops', {'-b', 'linewidth', 1}); hold on
shadedErrorBar(xvec, e_htom20, e_htom20_sem, ...
    'lineprops', {'-r', 'linewidth', 1});
ylim([8.5 13])
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('Wait time (s)');
title('\rm Expert 20: adapt to mixed')
ax2 = gca;
ax2.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off
axis square


% INFERENCE AGENT
nexttile
shadedErrorBar(xvec, inf_mtol, inf_mtol_sem, 'lineprops', ...
    {'-b', 'linewidth', 1})
hold on
shadedErrorBar(xvec, inf_mtoh, inf_mtoh_sem, ...
    'lineprops', {'-r', 'linewidth', 1});
ylim([9.5 13.5])
yl = ylim;
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('Wait time (s)');
title('\rm Inf: mixed to adapt')
ax5 = gca;
ax5.YRuler.TickLabelGapOffset = 1;
axis square
set(gca, 'TickDir', 'out'); box off

nexttile
shadedErrorBar(xvec, inf_ltom, inf_ltom_sem, ...
    'lineprops', {'-b', 'linewidth', 1}); hold on
shadedErrorBar(xvec, inf_htom, inf_htom_sem, ...
    'lineprops', {'-r', 'linewidth', 1});
ylim([9.5 13.5])
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('Wait time (s)');
title('\rm Inf: adapt to mixed')
ax6 = gca;
ax6.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off
axis square

% INFERENCE AGENT
nexttile
shadedErrorBar(xvec, inf_mtol20, inf_mtol20_sem, 'lineprops', ...
    {'-b', 'linewidth', 1})
hold on
shadedErrorBar(xvec, inf_mtoh20, inf_mtoh20_sem, ...
    'lineprops', {'-r', 'linewidth', 1});
ylim([9.5 13.5])
yl = ylim;
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('Wait time (s)');
title('\rm Inf 20: mixed to adapt')
ax5 = gca;
ax5.YRuler.TickLabelGapOffset = 1;
axis square
set(gca, 'TickDir', 'out'); box off

nexttile
shadedErrorBar(xvec, inf_ltom20, inf_ltom20_sem, ...
    'lineprops', {'-b', 'linewidth', 1}); hold on
shadedErrorBar(xvec, inf_htom20, inf_htom20_sem, ...
    'lineprops', {'-r', 'linewidth', 1});
ylim([9.5 13.5])
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('Wait time (s)');
title('\rm Inf 20: adapt to mixed')
ax6 = gca;
ax6.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off
axis square


set(gcf,'renderer','painter')
set(gcf, 'units', 'centimeters', 'position', fsize)

