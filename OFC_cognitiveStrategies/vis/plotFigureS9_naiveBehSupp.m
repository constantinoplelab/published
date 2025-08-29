function plotFigureS9_naiveBehSupp(processedBehaviorPath)
%Plot figure S9 - raw wait time block transition dynamics for the first 15
%   sessions and divisive normalization model predictions. Sigmoid fits to
%   z-scored 
% YOU MUST RUN processRawWaitTimes BEFORE RUNNING THIS FUNCTION

%INPUTS:
%   processedBehaviorPath = path to where you have saved the outputs from
%       running processRawWaitTimes

%load processed behavior data
load([processedBehaviorPath, filesep, 'naive.mat'])
load([processedBehaviorPath, filesep, 'divModel.mat'])

%general functions
sem = @(x) std(x,'omitnan')./sqrt(sum(any(~isnan(x), 2)));
setDefaultFigProps

twin = 40; %trial window for wait time dynamics plot
xvec = -twin:twin-1;
q = 1:4;
fsize = [10 10 20 12];

%% averages
n_mtol = mean(naive.mtolRaw, 'omitnan');
n_mtol_sem = sem(naive.mtolRaw);

n_mtoh = mean(naive.mtohRaw, 'omitnan');
n_mtoh_sem = sem(naive.mtohRaw);

n_ltom = mean(naive.ltomRaw, 'omitnan');
n_ltom_sem = sem(naive.ltomRaw);

n_htom = mean(naive.htomRaw, 'omitnan');
n_htom_sem = sem(naive.htomRaw);

n_mtol20 = mean(naive.mtol20Raw, 'omitnan');
n_mtol20_sem = sem(naive.mtol20Raw);

n_mtoh20 = mean(naive.mtoh20Raw, 'omitnan');
n_mtoh20_sem = sem(naive.mtoh20Raw);

n_ltom20 = mean(naive.ltom20Raw, 'omitnan');
n_ltom20_sem = sem(naive.ltom20Raw);

n_htom20 = mean(naive.htom20Raw, 'omitnan');
n_htom20_sem = sem(naive.htom20Raw);

div_mtol = mean(divModel.mtolRaw, 'omitnan');
div_mtol_sem = sem(divModel.mtolRaw);

div_mtoh = mean(divModel.mtohRaw, 'omitnan');
div_mtoh_sem = sem(divModel.mtohRaw);

div_ltom = mean(divModel.ltomRaw, 'omitnan');
div_ltom_sem = sem(divModel.ltomRaw);

div_htom = mean(divModel.htomRaw, 'omitnan');
div_htom_sem = sem(divModel.htomRaw);

div_mtol20 = mean(divModel.mtol20Raw, 'omitnan');
div_mtol20_sem = sem(divModel.mtol20Raw);

div_mtoh20 = mean(divModel.mtoh20Raw, 'omitnan');
div_mtoh20_sem = sem(divModel.mtoh20Raw);

div_ltom20 = mean(divModel.ltom20Raw, 'omitnan');
div_ltom20_sem = sem(divModel.ltom20Raw);

div_htom20 = mean(divModel.htom20Raw, 'omitnan');
div_htom20_sem = sem(divModel.htom20Raw);

%% fit sigmoids to individual rat transition dynamics
disp('Fitting sigmoid functions to transition dynamics')
disp('This will take ~4 minutes')

%window to fit sigmoid -- produces the most reasonable fits
wndw = [find(xvec==-10):find(xvec==25)];
nrats = length(naive.ltom);

finalParams_ltom_naive = nan(nrats, 4);
finalParams_mtol_naive = nan(nrats, 4);
finalParams_htom_naive = nan(nrats, 4);
finalParams_mtoh_naive = nan(nrats, 4);
for rr = 1:nrats
    try
        finalParams_ltom_naive(rr, :) = my_fit_sigmoid(xvec(wndw), naive.ltom(rr, wndw), 50);
        finalParams_mtol_naive(rr, :) = my_fit_sigmoid(xvec(wndw), naive.mtol(rr, wndw), 50);
        finalParams_htom_naive(rr, :) = my_fit_sigmoid(xvec(wndw), naive.htom(rr, wndw), 50);
        finalParams_mtoh_naive(rr, :) = my_fit_sigmoid(xvec(wndw), naive.mtoh(rr, wndw), 50);

    catch
        continue
    end
end

% Fit sigmoid to naive data -- aligned to incongruent

finalParams_ltomIncong_naive = nan(nrats, 4);
finalParams_htomIncong_naive = nan(nrats, 4);
for rr = 1:nrats
    try
        finalParams_ltomIncong_naive(rr, :) = my_fit_sigmoid(xvec(wndw), ...
            naive.ltom_incong(rr, wndw), 50);  
        finalParams_htomIncong_naive(rr, :) = my_fit_sigmoid(xvec(wndw), ...
            naive.htom_incong(rr, wndw), 50);

    catch
        continue
    end
end

%% sigmoid fit averages and stats
n_mtol_fit = median(finalParams_mtol_naive(:,3), 'omitnan');
n_ltom_fit = median(finalParams_ltom_naive(:,3), 'omitnan');
n_mtoh_fit = median(finalParams_mtoh_naive(:,3), 'omitnan'); 
n_htom_fit = median(finalParams_htom_naive(:,3), 'omitnan');
n_ltomIncong_fit = median(finalParams_ltomIncong_naive(:,3), 'omitnan');
n_htomIncong_fit = median(finalParams_htomIncong_naive(:,3), 'omitnan');

p_low_n = signrank(finalParams_mtol_naive(:,3), finalParams_ltom_naive(:,3));
p_high_n = signrank(finalParams_mtoh_naive(:,3), finalParams_htom_naive(:,3));

%% Plot First 15 sessions and divisive norm. model

figure; hold on
tiledlayout(3, 4, 'TileSpacing', 'compact')

%FIRST 15 SESSIONS, congruent volumes
nexttile
shadedErrorBar(xvec, n_mtol, n_mtol_sem, ...
    'lineprops', {'b', 'linewidth', 1});
shadedErrorBar(xvec, n_mtoh, n_mtoh_sem, ...
    'lineprops', {'r', 'linewidth', 1});
ylim([9.5 13.5])
yl = ylim;
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('Wait time (s)');
title('\rm Naive: mixed to adapt')
ax1 = gca;
ax1.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off
axis square

nexttile
shadedErrorBar(xvec, n_ltom, n_ltom_sem, ...
    'lineprops', {'b', 'linewidth', 1});
shadedErrorBar(xvec, n_htom, n_htom_sem, ...
    'lineprops', {'r', 'linewidth', 1});
ylim([9.5 13.5])
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('Wait time (s)');
title('\rm Naive: adapt to mixed')
ax2 = gca;
ax2.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off
axis square

%FIRST 15 SESSIONS
nexttile
shadedErrorBar(xvec, n_mtol20, n_mtol20_sem, ...
    'lineprops', {'b', 'linewidth', 1});
shadedErrorBar(xvec, n_mtoh20, n_mtoh20_sem, ...
    'lineprops', {'r', 'linewidth', 1});
ylim([9.5 13.5])
yl = ylim;
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('Wait time (s)');
title('\rm Naive 20: mixed to adapt')
ax3 = gca;
ax3.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off
axis square

nexttile
shadedErrorBar(xvec, n_ltom20, n_ltom20_sem, ...
    'lineprops', {'b', 'linewidth', 1});
shadedErrorBar(xvec, n_htom20, n_htom20_sem, ...
    'lineprops', {'r', 'linewidth', 1});
ylim([9.5 13.5])
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('Wait time (s)');
title('\rm Naive 20: adapt to mixed')
ax4 = gca;
ax4.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off
axis square


nexttile
shadedErrorBar(xvec, div_mtol, div_mtol_sem, ...
    'lineprops', {'b', 'linewidth', 1});
shadedErrorBar(xvec, div_mtoh, div_mtoh_sem, ...
    'lineprops', {'r', 'linewidth', 1});
ylim([2 9])
yl = ylim;
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('Wait time (s)');
title('\rm Div: mixed to adapt')
ax7 = gca;
ax7.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off
axis square

%DIVISIVE NORM MODEL, congruent volumes
nexttile
shadedErrorBar(xvec, div_ltom, div_ltom_sem, ...
    'lineprops', {'b', 'linewidth', 1});
shadedErrorBar(xvec, div_htom, div_htom_sem, ...
    'lineprops', {'r', 'linewidth', 1});
ylim([2 9])
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('Wait time (s)');
title('\rm Div: adapt to mixed')
ax8 = gca;
ax8.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off
axis square

nexttile
shadedErrorBar(xvec, div_mtol20, div_mtol20_sem, ...
    'lineprops', {'b', 'linewidth', 1});
shadedErrorBar(xvec, div_mtoh20, div_mtoh20_sem, ...
    'lineprops', {'r', 'linewidth', 1});
ylim([2 9])
yl = ylim;
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('Wait time (s)');
title('\rm Div 20: mixed to adapt')
ax7 = gca;
ax7.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off
axis square

%DIVISIVE NORM MODEL, 20ul only
nexttile
shadedErrorBar(xvec, div_ltom20, div_ltom20_sem, ...
    'lineprops', {'b', 'linewidth', 1});
shadedErrorBar(xvec, div_htom20, div_htom20_sem, ...
    'lineprops', {'r', 'linewidth', 1});
ylim([2 9])
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('Wait time (s)');
title('\rm Div 20: adapt to mixed')
ax8 = gca;
ax8.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off
axis square

% plot sigmoid fits for mixed to low and low to mixed transitions
nexttile
histogram(finalParams_mtol_naive(:,3), 'facecolor', 'b')
yl = ylim;
hold on
xline(n_mtol_fit, '--b')
histogram(finalParams_ltom_naive(:,3), 'facecolor', 'k')
xline(n_ltom_fit, '--k')
set(gca, 'TickDir', 'out'); box off;
xlabel('X offset')
ylabel('# rats')
title(['\rm p = ' string(p_low_n)])
axis square

% plot sigmoid fits for mixed to high and high to mixed transitions
nexttile
histogram(finalParams_mtoh_naive(:,3), 'facecolor', 'r')
yl = ylim;
hold on
xline(n_mtoh_fit, '--r')
histogram(finalParams_htom_naive(:,3), 'facecolor', 'k')
xline(n_htom_fit, '--k')
set(gca, 'TickDir', 'out'); box off;
xlabel('X offset')
ylabel('# rats')
title(['\rm p = ' string(p_high_n)])
axis square

% plot sigmoid fits for low to mixed transitions, aligned to the first
% incongruent trial
nexttile
hold on
histogram(finalParams_ltomIncong_naive(:,3), 'facecolor', 'k')
xline(n_ltomIncong_fit, '--k')
set(gca, 'TickDir', 'out'); box off;
xlabel('X offset')
ylabel('# rats')
title(['\rm' string(n_ltomIncong_fit)])
axis square

% plot sigmoid fits for high to mixed transitions, aligned to the first
% incongruent trial
nexttile
hold on
histogram(finalParams_htomIncong_naive(:,3), 'facecolor', 'k')
xline(n_htomIncong_fit, '--k')
set(gca, 'TickDir', 'out'); box off;
xlabel('X offset')
ylabel('# rats')
title(['\rm' string(n_htomIncong_fit)])
axis square


set(gcf,'renderer','painter')
set(gcf, 'units', 'centimeters', 'position', fsize)

