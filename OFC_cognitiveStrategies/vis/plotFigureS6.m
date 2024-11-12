function plotFigureS6(behaviorPath)
%plots data for figure S6

%INPUTS:
%   behaviorPath = path to behavior data from all 349 rats

%paths
load(strcat(behaviorPath, 'ratList.mat')); %list of rat data to use for model simulations

A = load(strcat(behaviorPath, 'ratTrial_S026.mat'));
ratTrial = A.A;
ratTrial.tau = 2.5;

twin = 40; %trial window to use for wait time dynamics simulations
binSize = 5; %number of trials to bin over
smoothfactor = 5;

x = 1:5;
xvec = -twin:5:twin-5; %for 5 trial window, first bin after block change is 0

%general functions
sem = @(x) std(x,'omitnan')./sqrt(sum(any(~isnan(x), 2)));
setDefaultFigProps

%% LAMBDA SIMULATIONS

%simulate wait times for an example rat -- manipulate the quality of the
%prior
[simC_lambda, simM_lambda] = simWTcurves(ratTrial, 'lambda');

%simulate transition dynamics averaged over all rats
[simControl_lambda, simMuscimol_lambda] = simDynamics(ratList, behaviorPath, ...
    'lambda', twin, binSize, smoothfactor);


lambda_mtol = mean(simMuscimol_lambda.mtol, 'omitnan');
lambda_mtol_sem = sem(simMuscimol_lambda.mtol);

lambda_mtoh = mean(simMuscimol_lambda.mtoh, 'omitnan');
lambda_mtoh_sem = sem(simMuscimol_lambda.mtoh);

opt_mtol = mean(simControl_lambda.mtol, 'omitnan');
opt_mtol_sem = sem(simControl_lambda.mtol);

opt_mtoh = mean(simControl_lambda.mtoh, 'omitnan');
opt_mtoh_sem = sem(simControl_lambda.mtoh);

p_optEarly = signrank(simControl_lambda.mtol(:,12), simControl_lambda.mtoh(:,12));
p_optLate = signrank(simControl_lambda.mtol(:,16), simControl_lambda.mtoh(:,16));
p_lambdaEarly = signrank(simMuscimol_lambda.mtol(:,12), simMuscimol_lambda.mtoh(:,12));
p_lambdaLate = signrank(simMuscimol_lambda.mtol(:,16), simMuscimol_lambda.mtoh(:,16));

%% KAPPA SIMULATIONS

%simulate wait times for an example rat -- manipulate how similar the
%opportunity costs are for each block
[simC_kappa, simM_kappa] = simWTcurves(ratTrial, 'kappa');

%simulate transition dynamics averaged over all rats
[simControl_kappa, simMuscimol_kappa] = simDynamics(ratList, behaviorPath, ...
    'kappa', twin, binSize, smoothfactor);


kappa_mtol = mean(simMuscimol_kappa.mtol, 'omitnan');
kappa_mtol_sem = sem(simMuscimol_kappa.mtol);

kappa_mtoh = mean(simMuscimol_kappa.mtoh, 'omitnan');
kappa_mtoh_sem = sem(simMuscimol_kappa.mtoh);

optK_mtol = mean(simControl_kappa.mtol, 'omitnan');
optK_mtol_sem = sem(simControl_kappa.mtol);

optK_mtoh = mean(simControl_kappa.mtoh, 'omitnan');
optK_mtoh_sem = sem(simControl_kappa.mtoh);

p_optKEarly = signrank(simControl_kappa.mtol(:,12), simControl_kappa.mtoh(:,12));
p_optKLate = signrank(simControl_kappa.mtol(:,16), simControl_kappa.mtoh(:,16));
p_kappaEarly = signrank(simMuscimol_kappa.mtol(:,12), simMuscimol_kappa.mtoh(:,12));
p_kappaLate = signrank(simMuscimol_kappa.mtol(:,16), simMuscimol_kappa.mtoh(:,16));

%% model schematics

tau = 2.5;
C = 1 - 0.15;
t = 0:0.5:20;
%model value function
v = @(tau, t, C) (1/tau)*((C*exp(-t/tau)) / (1 - C + C*exp(-t/tau)));
value = arrayfun(@(x) v(tau, t(x), C), 1:length(t));

% opportunity cost sim
rng(23)

r = [randsample(2:6, 40, true),...
    randsample(2:4, 40, true),...
    randsample(2:6, 40, true),...
    randsample(4:6, 40, true),...
    randsample(2:6, 40, true)]';
r(end-20:end-15) = randsample(2:4, 6, true);

rvec = [5 10 20 40 80];
A_model.reward = rvec(r-1)';
A_model.wait_time = ones(size(A_model.reward));
A_model.prob_catch = zeros(size(r));
A_model.tau = tau;
A_model.ntrials = length(r);

params = [0.23 0.28 0.18 0.13];
[~, ~, ~, Belief, ~, Kappa] = GenerateSynthData_Bayes(params, A_model, ...
    'logn', 1, 8);

%%

figure;
tiledlayout(3, 4, 'TileSpacing', 'tight')

%model value function
nexttile
plot(t, value, 'k', 'linewidth', 1.5)
yline(0.2, '--r', 'linewidth', 1)
yline(0.08, '--b', 'linewidth', 1)
ylim([0 0.40])
yticks([])
xticks([])
box off
xlabel('Time')
ylabel('Offer value')
title('\rm Trial value function')
ax1 = gca;
ax1.YRuler.TickLabelGapOffset = 1;
axis square

% example model predicted wait time curve
nexttile
plot(x, simC_lambda.hi, '-r', 'linewidth', 1.5); hold on
plot(x, simC_lambda.lo, '-b', 'linewidth', 1.5)
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
ylim([10 14])
xlabel('Reward offer')
ylabel('Wait wime (s)')
ax2 = gca;
ax2.YRuler.TickLabelGapOffset = 1;
axis square
title('\rm Model')

% example model predicted posterior beliefs during a session
nexttile
yl = [0 1];
fill([0 40 40 0], [yl(1) yl(1) yl(2) yl(2)],...
    'k', 'facealpha', 0.15, 'edgecolor', 'none'); hold on
fill([40 80 80 40], [yl(1) yl(1) yl(2) yl(2)],...
    'b', 'facealpha', 0.15, 'edgecolor', 'none')
fill([80 120 120 80], [yl(1) yl(1) yl(2) yl(2)],...
    'k', 'facealpha', 0.15, 'edgecolor', 'none')
fill([120 160 160 120], [yl(1) yl(1) yl(2) yl(2)],...
    'r', 'facealpha', 0.15, 'edgecolor', 'none')
fill([160 200 200 160], [yl(1) yl(1) yl(2) yl(2)],...
    'k', 'facealpha', 0.15, 'edgecolor', 'none')

plot(Belief(1,:), 'k', 'linewidth', 0.75)
plot(Belief(2,:), 'r', 'linewidth', 0.75)
plot(Belief(3,:), 'b', 'linewidth', 0.75)
ylim(yl)
yticks([])
xticks([])
box off
xlabel('Trial number')
ylabel('Posterior belief')
ax3 = gca;
ax3.YRuler.TickLabelGapOffset = 1;
axis square

% example model predicted opportunity costs during a session
% the opportunity cost is selected on each trial based on the maximum 
% posterior belief shown in panel c

nexttile
yl = [0.15 0.28];
fill([0 40 40 0], [yl(1) yl(1) yl(2) yl(2)],...
    'k', 'facealpha', 0.15, 'edgecolor', 'none'); hold on
fill([40 80 80 40], [yl(1) yl(1) yl(2) yl(2)],...
    'b', 'facealpha', 0.15, 'edgecolor', 'none')
fill([80 120 120 80], [yl(1) yl(1) yl(2) yl(2)],...
    'k', 'facealpha', 0.15, 'edgecolor', 'none')
fill([120 160 160 120], [yl(1) yl(1) yl(2) yl(2)],...
    'r', 'facealpha', 0.15, 'edgecolor', 'none')
fill([160 200 200 160], [yl(1) yl(1) yl(2) yl(2)],...
    'k', 'facealpha', 0.15, 'edgecolor', 'none')

plot(Kappa, '-k', 'linewidth', 1)
ylim(yl)
yticks([])
xticks([])
box off
xlabel('Trial number')
ylabel('Opportunity cost')
ax4 = gca;
ax4.YRuler.TickLabelGapOffset = 1;
axis square

%manipulate prior quality
nexttile
plot(x, simC_lambda.hi, '-r', 'linewidth', 1.5); hold on
plot(x, simC_lambda.lo, '-b', 'linewidth', 1.5)
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
ylim([10 14])
xlabel('Reward offer')
ylabel('Mean wait wime (s)')
ax5 = gca;
ax5.YRuler.TickLabelGapOffset = 1;
axis square
title('\rm Quality prior')

nexttile
plot(x, simM_lambda.hi, 'color', [1 0.6 0.6], 'linewidth', 1.5); hold on
plot(x, simM_lambda.lo, 'color', [0.6 0.6 1], 'linewidth', 1.5)
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
ylim([10 14])
xlabel('Reward offer')
ylabel('Mean wait wime (s)')
ax6 = gca;
ax6.YRuler.TickLabelGapOffset = 1;
axis square
linkaxes([ax5 ax6])
title('\rm Poor prior')

nexttile
shadedErrorBar(xvec, opt_mtol, opt_mtol_sem, ...
    'lineprops', {'-b', 'linewidth', 1.5});
shadedErrorBar(xvec, opt_mtoh, opt_mtoh_sem, ...
    'lineprops', {'-r', 'linewidth', 1.5});
shadedErrorBar(xvec, lambda_mtol, lambda_mtol_sem, ...
    'lineprops', {'color', [0.6 0.6 1], 'linewidth', 1.5});
shadedErrorBar(xvec, lambda_mtoh, lambda_mtoh_sem, ...
    'lineprops', {'color', [1 0.6 0.6], 'linewidth', 1.5});
yl = [-2 2];
xlim([-30 35]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('\Delta z-scored wait time');
title('\rm Mixed into adaptation');
axis square
ax7 = gca;
ax7.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off

nexttile
errorbar([0.25 0.5], [opt_mtol(12) opt_mtol(16)], ...
    [opt_mtol_sem(12) opt_mtol_sem(16)], '_b', 'markersize', 15, ...
    'linewidth', 1, 'linestyle', 'none'); hold on
errorbar([0.25 0.5], [opt_mtoh(12) opt_mtoh(16)], ...
    [opt_mtoh_sem(12) opt_mtoh_sem(16)], '_r', 'markersize', 15, ...
    'linewidth', 1, 'linestyle', 'none'); 
errorbar([1 1.25], [lambda_mtol(12) lambda_mtol(16)], ...
    [lambda_mtol_sem(12) lambda_mtol_sem(16)], 'marker', '_', 'color', ...
    [0.6 0.6 1], 'markersize', 15, 'linewidth', 1, 'linestyle', 'none');
errorbar([1 1.25], [lambda_mtoh(12) lambda_mtoh(16)], ...
    [lambda_mtoh_sem(12) lambda_mtoh_sem(16)], 'marker', '_', 'color', ...
    [1 0.6 0.6], 'markersize', 15, 'linewidth', 1, 'linestyle', 'none'); hold on
yline(0, '--k', 'linewidth', 1)
xlim([0 1.5])
ylim([-1.5 1.5])
ylabel('\Delta z-scored wait time');
% title({strcat(string(p_optEarly), {'  '}, string(p_lambdaEarly)), ...
%     strcat(string(p_optLate), {'  '}, string(p_lambdaLate))})
axis square
ax8 = gca;
ax8.YRuler.TickLabelGapOffset = 1;
set(gca, 'xTick', [0.25 0.5 1 1.25]);
set(gca, 'XTickLabels', {'Early' 'Late' 'Early' 'Late'});
set(gca, 'TickDir', 'out'); box off   

%manipulate opportunity costs
nexttile
plot(x, simC_kappa.hi, '-r', 'linewidth', 1.5); hold on
plot(x, simC_kappa.lo, '-b', 'linewidth', 1.5)
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
ylim([10 14])
xlabel('Reward offer')
ylabel('Mean wait time (s)')
ax9 = gca;
ax9.YRuler.TickLabelGapOffset = 1;
axis square
title('\rm Distinct kappas')

nexttile
plot(x, simM_kappa.hi, 'color', [1 0.6 0.6], 'linewidth', 1.5); hold on
plot(x, simM_kappa.lo, 'color', [0.6 0.6 1], 'linewidth', 1.5)
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
ylim([10 14])
xlabel('Reward offer')
ylabel('Mean wait time (s)')
ax10 = gca;
ax10.YRuler.TickLabelGapOffset = 1;
axis square
linkaxes([ax9 ax10])
title('\rm Similar kappas')

nexttile
shadedErrorBar(xvec, optK_mtol, optK_mtol_sem, ...
    'lineprops', {'-b', 'linewidth', 1.5});
shadedErrorBar(xvec, optK_mtoh, optK_mtoh_sem, ...
    'lineprops', {'-r', 'linewidth', 1.5});
shadedErrorBar(xvec, kappa_mtol, kappa_mtol_sem, ...
    'lineprops', {'color', [0.6 0.6 1], 'linewidth', 1.5});
shadedErrorBar(xvec, kappa_mtoh, kappa_mtoh_sem, ...
    'lineprops', {'color', [1 0.6 0.6], 'linewidth', 1.5});
yl = [-2 2];
xlim([-30 35]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('\Delta z-scored wait time');
title('\rm Mixed into adaptation');
axis square
ax11 = gca;
ax11.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off

nexttile
errorbar([0.25 0.5], [optK_mtol(12) optK_mtol(16)], ...
    [optK_mtol_sem(12) optK_mtol_sem(16)], '_b', 'markersize', 15, ...
    'linewidth', 1, 'linestyle', 'none'); hold on
errorbar([0.25 0.5], [optK_mtoh(12) optK_mtoh(16)], ...
    [optK_mtoh_sem(12) optK_mtoh_sem(16)], '_r', 'markersize', 15, ...
    'linewidth', 1, 'linestyle', 'none'); 
errorbar([1 1.25], [kappa_mtol(12) kappa_mtol(16)], ...
    [kappa_mtol_sem(12) kappa_mtol_sem(16)], 'marker', '_', 'color', ...
    [0.6 0.6 1], 'markersize', 15, 'linewidth', 1, 'linestyle', 'none');
errorbar([1 1.25], [kappa_mtoh(12) kappa_mtoh(16)], ...
    [kappa_mtoh_sem(12) kappa_mtoh_sem(16)], 'marker', '_', 'color', ...
    [1 0.6 0.6], 'markersize', 15, 'linewidth', 1, 'linestyle', 'none'); hold on
yline(0, '--k', 'linewidth', 1)
xlim([0 1.5])
ylim([-1.5 1.5])
ylabel('\Delta z-scored wait time');
% title({strcat(string(p_optKEarly), {'  '}, string(p_kappaEarly)), ...
%     strcat(string(p_optKLate), {'  '}, string(p_kappaLate))})
axis square
ax12 = gca;
ax12.YRuler.TickLabelGapOffset = 1;
set(gca, 'xTick', [0.25 0.5 1 1.25]);
set(gca, 'XTickLabels', {'Early' 'Late' 'Early' 'Late'});
set(gca, 'TickDir', 'out'); box off   
set(gcf,'renderer','painter')
set(gcf, 'units', 'centimeters', 'position', [10 10 20 14])

