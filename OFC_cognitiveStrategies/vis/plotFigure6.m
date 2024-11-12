function plotFigure6(muscimolPath, behaviorPath)

%INPUTS: 
%   muscimolPath = path to downloaded muscimol data. Folder contains the
%       subfolders BehaviorData and Physiology
%   behaviorPath = path to behavior data from all 349 rats

muscimolBehPath = [muscimolPath, filesep, 'BehaviorData', filesep];


ratListC = {'ratTrial_S001Control'; 'ratTrial_S002Control'; ...
    'ratTrial_S003Control'; 'ratTrial_S009Control'; 'ratTrial_S014Control'; ...
    'ratTrial_S017Control'; 'ratTrial_S018Control'; 'ratTrial_S027Control'; ...
    'ratTrial_S028Control'};

ratListM = {'ratTrial_S001Muscimol'; 'ratTrial_S002Muscimol'; ...
    'ratTrial_S003Muscimol'; 'ratTrial_S009Muscimol'; ...
    'ratTrial_S014Muscimol'; 'ratTrial_S017Muscimol'; 'ratTrial_S018Muscimol'; ...
    'ratTrial_S027Muscimol'; 'ratTrial_S028Muscimol'};

load(strcat(behaviorPath, 'ratList.mat')); %list of rat data to use for model simulations

%load example rat data for model wait time simulations
A = load(strcat(behaviorPath, 'ratTrial_S026.mat'));
ratTrial = A.A;
ratTrial.tau = 2.5; %time constant of the exponential delay distribution

%inputs for transition dynamics computations
twin = 40; %window used for block transition dynamics, adaptation blocks are 40 successful trials
binSize = 5; %average block transiont dynamics over bins of 5 trials
smoothfactor = 5;

%general functions
sem = @(x) std(x,'omitnan')./sqrt(sum(any(~isnan(x), 2)));
setDefaultFigProps
fsize = [10 10 20 13];

x = 1:5;
xvec = -twin:5:twin-5;

%% Process muscimol rat data
[control, muscimol] = processMuscimolData(ratListC, ratListM, ...
    muscimolBehPath, twin, binSize, smoothfactor);

%%
avgCh = mean(control.hi);
semCh = sem(control.hi);

avgCl = mean(control.lo);
semCl = sem(control.lo);

avgMh = mean(muscimol.hi);
semMh = sem(muscimol.hi);

avgMl = mean(muscimol.lo);
semMl = sem(muscimol.lo);

c_mtol = mean(control.mtol, 'omitnan');
c_mtol_sem = sem(control.mtol);

c_mtoh = mean(control.mtoh, 'omitnan');
c_mtoh_sem = sem(control.mtoh);

m_mtol = mean(muscimol.mtol, 'omitnan');
m_mtol_sem = sem(muscimol.mtol);

m_mtoh = mean(muscimol.mtoh, 'omitnan');
m_mtoh_sem = sem(muscimol.mtoh);

p_cEarly = signrank(control.mtol(:,12), control.mtoh(:,12));
p_cLate = signrank(control.mtol(:,16), control.mtoh(:,16));
p_mEarly = signrank(muscimol.mtol(:,12), muscimol.mtoh(:,12));
p_mLate = signrank(muscimol.mtol(:,16), muscimol.mtoh(:,16));

%% LAPSE RATE SIMULATIONS
[simC, simM] = simWTcurves(ratTrial, 'lapse');

[simControl_lapse, simMuscimol_lapse] = simDynamics(ratList, ...
    behaviorPath, 'lapse', twin, binSize, smoothfactor);
 
%%

lapse_mtol = mean(simMuscimol_lapse.mtol, 'omitnan');
lapse_mtol_sem = sem(simMuscimol_lapse.mtol);

lapse_mtoh = mean(simMuscimol_lapse.mtoh, 'omitnan');
lapse_mtoh_sem = sem(simMuscimol_lapse.mtoh);

opt_mtol = mean(simControl_lapse.mtol, 'omitnan');
opt_mtol_sem = sem(simControl_lapse.mtol);

opt_mtoh = mean(simControl_lapse.mtoh, 'omitnan');
opt_mtoh_sem = sem(simControl_lapse.mtoh);

p_optEarly = signrank(simControl_lapse.mtol(:,12), simControl_lapse.mtoh(:,12));
p_optLate = signrank(simControl_lapse.mtol(:,16), simControl_lapse.mtoh(:,16));
p_lapseEarly = signrank(simMuscimol_lapse.mtol(:,12), simMuscimol_lapse.mtoh(:,12));
p_lapseLate = signrank(simMuscimol_lapse.mtol(:,16), simMuscimol_lapse.mtoh(:,16));

%% plot

figure; hold on
tiledlayout(2, 4, 'TileSpacing', 'compact')

%wait time curves
nexttile
shadedErrorBar(x, avgCh, semCh, 'lineprops', {'-r', 'linewidth', 1.5})
shadedErrorBar(x, avgCl, semCl, 'lineprops', {'-b', 'linewidth', 1.5})
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
xlabel('Reward offer')
ylabel('Mean wait wime (s)')
ax1 = gca;
ax1.YRuler.TickLabelGapOffset = 1;
axis square
title('\rm Control')

nexttile
shadedErrorBar(x, avgMh, semMh, 'lineprops', {'color', [1 0.6 0.6], ...
    'linewidth', 1.5})
shadedErrorBar(x, avgMl, semMl, 'lineprops', {'color', [0.6 0.6 1],...
    'linewidth', 1.5})
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
xlabel('Reward offer')
ylabel('Mean wait time (s)')
ax2 = gca;
ax2.YRuler.TickLabelGapOffset = 1;
axis square
title('\rm Muscimol')
set(gcf,'renderer','painter')
set(gcf, 'Color', [1 1 1]);

%wait time dynamics - muscimol
nexttile
shadedErrorBar(xvec, c_mtol, c_mtol_sem, ...
    'lineprops', {'-b', 'linewidth', 1.5}); hold on
shadedErrorBar(xvec, c_mtoh, c_mtoh_sem, ...
    'lineprops', {'-r', 'linewidth', 1.5});
shadedErrorBar(xvec, m_mtol, m_mtol_sem, ...
    'lineprops', {'color', [0.6 0.6 1], 'linewidth', 1.5});
shadedErrorBar(xvec, m_mtoh, m_mtoh_sem, ...
    'lineprops', {'color', [1 0.6 0.6], 'linewidth', 1.5});
ylim([-0.4 0.4])
yl = ylim;
xlim([-30 35]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('\Delta z-scored wait time');
title('\rm Mixed into adaptation');
axis square
ax3 = gca;
ax3.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off

%early and late in a block
nexttile
errorbar([0.25 0.5], [c_mtol(12) c_mtol(16)], ...
    [c_mtol_sem(12) c_mtol_sem(16)], '_b', 'markersize', 15, 'linewidth', 1, ...
    'linestyle', 'none'); hold on
errorbar([0.25 0.5], [c_mtoh(12) c_mtoh(16)], ...
    [c_mtoh_sem(12) c_mtoh_sem(16)], '_r', 'markersize', 15, 'linewidth', 1, ...
    'linestyle', 'none'); 
errorbar([1 1.25], [m_mtol(12) m_mtol(16)], ...
    [m_mtol_sem(12) m_mtol_sem(16)], 'marker', '_', 'color', [0.6 0.6 1], ...
    'markersize', 15, 'linewidth', 1, 'linestyle', 'none');
errorbar([1 1.25], [m_mtoh(12) m_mtoh(16)], ...
    [m_mtoh_sem(12) m_mtoh_sem(16)], 'marker', '_', 'color', [1 0.6 0.6], ...
    'markersize', 15, 'linewidth', 1, 'linestyle', 'none'); hold on
yline(0, '--k', 'linewidth', 1)
xlim([0 1.5])
ylim([-0.4 0.4])
ylabel('\Delta z-scored wait time');
title({strcat(string(p_cEarly), {'  '}, string(p_mEarly)), ...
    strcat(string(p_cLate), {'  '}, string(p_mLate))})
axis square
ax4 = gca;
ax4.YRuler.TickLabelGapOffset = 1;
set(gca, 'xTick', [0.25 0.5 1 1.25]);
set(gca, 'XTickLabels', {'Early' 'Late' 'Early' 'Late'});
set(gca, 'TickDir', 'out'); box off   
set(gcf,'renderer','painter')
set(gcf, 'units', 'centimeters', 'position', fsize)

% wait times - model
nexttile
plot(x, simC.hi, '-r', 'linewidth', 1.5); hold on
plot(x, simC.lo, '-b', 'linewidth', 1.5)
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
xlabel('Reward offer')
ylabel('Mean wait time (s)')
ax5 = gca;
ax5.YRuler.TickLabelGapOffset = 1;
axis square
title('\rm Low lapse')

nexttile
plot(x, simM.hi, 'color', [1 0.6 0.6], 'linewidth', 1.5); hold on
plot(x, simM.lo, 'color', [0.6 0.6 1], 'linewidth', 1.5)
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
xlabel('Reward offer')
ylabel('Mean wait time (s)')
ax6 = gca;
ax6.YRuler.TickLabelGapOffset = 1;
axis square
linkaxes([ax1 ax2 ax5 ax6])
title('\rm High lapse')

%dynamics - model
nexttile
shadedErrorBar(xvec, opt_mtol, opt_mtol_sem, ...
    'lineprops', {'-b', 'linewidth', 1.5});
shadedErrorBar(xvec, opt_mtoh, opt_mtoh_sem, ...
    'lineprops', {'-r', 'linewidth', 1.5});
shadedErrorBar(xvec, lapse_mtol, lapse_mtol_sem, ...
    'lineprops', {'color', [0.6 0.6 1], 'linewidth', 1.5});
shadedErrorBar(xvec, lapse_mtoh, lapse_mtoh_sem, ...
    'lineprops', {'color', [1 0.6 0.6], 'linewidth', 1.5});
ylim([-1.5 1.5])
yl = ylim;
xlim([-30 35]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('\Delta z-scored wait time');
title('\rm Mixed into adaptation');
axis square
ax7 = gca;
ax7.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off

%early and late - model
nexttile
errorbar([0.25 0.5], [opt_mtol(12) opt_mtol(16)], ...
    [opt_mtol_sem(12) opt_mtol_sem(16)], '_b', 'markersize', 15, ...
    'linewidth', 1, 'linestyle', 'none'); hold on
errorbar([0.25 0.5], [opt_mtoh(12) opt_mtoh(16)], ...
    [opt_mtoh_sem(12) opt_mtoh_sem(16)], '_r', 'markersize', 15, ...
    'linewidth', 1, 'linestyle', 'none'); 
errorbar([1 1.25], [lapse_mtol(12) lapse_mtol(16)], ...
    [lapse_mtol_sem(12) lapse_mtol_sem(16)], 'marker', '_', 'color', ...
    [0.6 0.6 1], 'markersize', 15, 'linewidth', 1, 'linestyle', 'none');
errorbar([1 1.25], [lapse_mtoh(12) lapse_mtoh(16)], ...
    [lapse_mtoh_sem(12) lapse_mtoh_sem(16)], 'marker', '_', 'color', ...
    [1 0.6 0.6], 'markersize', 15, 'linewidth', 1, 'linestyle', 'none'); hold on
yline(0, '--k', 'linewidth', 1)
xlim([0 1.5])
ylim([-1.5 1.5])
ylabel('\Delta z-scored wait time');
title({strcat(string(p_optEarly), {'  '}, string(p_lapseEarly)), ...
    strcat(string(p_optLate), {'  '}, string(p_lapseLate))})
axis square
ax8 = gca;
ax8.YRuler.TickLabelGapOffset = 1;
set(gca, 'xTick', [0.25 0.5 1 1.25]);
set(gca, 'XTickLabels', {'Early' 'Late' 'Early' 'Late'});
set(gca, 'TickDir', 'out'); box off   
set(gcf,'renderer','painter')
set(gcf, 'units', 'centimeters', 'position', fsize)

