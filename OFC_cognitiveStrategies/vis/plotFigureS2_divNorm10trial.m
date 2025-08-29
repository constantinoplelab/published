function plotFigureS2_divNorm10trial(behaviorPath)
% Process and plot data for figure S2. Plots divisive normalization model 
% predictions for a trial integration window of 10. Plots conditional wait 
% times for 20ul trials for experts, inference model, and 10 trial divisive
% normalization model. Plot sigmoid fit histograms for 10 trial divisive 
% normalization model and inference model

%INPUTS: 
%   behaviorPath = path to behavior data from the repository

R = load(strcat(behaviorPath, 'ratList.mat')); %list of rat data to use for behavior plots and model simulations
ratList = R.ratList;
nrats = length(ratList);

behaviorPath_expert = behaviorPath;

sumtrials = 10; %number of trials to sum over for previous sum analysis

twin = 40; %trial window for wait time dynamics plot
xvec = -twin:twin-1;
smoothfactor = 5;
binSize = 1;
q = 1:4;

sem = @(x) std(x,'omitnan')./sqrt(sum(any(~isnan(x), 2)));
setDefaultFigProps

%%
disp('Running simulations and processing expert rat data')
disp('This will take ~ 3 minutes')

for rr = 1:nrats
    
    E = load(strcat(behaviorPath_expert, 'ratTrial_', ratList{rr}, '.mat'));

    % simulate divisive normalization agent with smaller integration window
    div = E.A;
        % linear rewards - run first so all other analyses use wait times
        % based on log2 rewards as in the rest of the manuscript
    divnorm_wt_linear = divisiveNorm_fun_2(div, [50, 0.15, sumtrials], 'linear');
    div.wait_time = divnorm_wt_linear;
    [divModel.ltom_linear(rr,:), divModel.htom_linear(rr,:), divModel.mtol_linear(rr,:), ...
        divModel.mtoh_linear(rr,:), ~, ~, ~] =...
        block_dynamics_wt_binTrials(div, twin, binSize, smoothfactor);
        %  log2 rewards
    divnorm_wt = divisiveNorm_fun_2(div, [50, 0.15, sumtrials], 'log');
    div.wait_time = divnorm_wt;
    [divModel.ltom(rr,:), divModel.htom(rr,:), divModel.mtol(rr,:), ...
        divModel.mtoh(rr,:), ~, ~, ~] =...
        block_dynamics_wt_binTrials(div, twin, binSize, smoothfactor);


    % simulate inference agent
    inf = E.A;
        % linear rewards
    [inf_wt_linear, ~, ~] = GenerateSynthData_Bayes_SS([0.25 0.3 0.2 .13], ...
        inf, 'logn', 1, 8, 'linear');
    inf.wait_time = inf_wt_linear;
    [infModel.ltom_linear(rr,:), infModel.htom_linear(rr,:), infModel.mtol_linear(rr,:), ...
        infModel.mtoh_linear(rr,:), ~, ~, ~] =...
        block_dynamics_wt_binTrials(inf, twin, binSize, smoothfactor);    
        % log2 rewards
    [inf_wt, inf_wt_mdl, BlkInf] = GenerateSynthData_Bayes_SS([0.25 0.3 0.2 .13], ...
        inf, 'logn', 1, 8, 'log');
    inf.wait_time = inf_wt; 
    [infModel.ltom(rr,:), infModel.htom(rr,:), infModel.mtol(rr,:), ...
        infModel.mtoh(rr,:), ~, ~, ~] =...
        block_dynamics_wt_binTrials(inf, twin, binSize, smoothfactor);


    % split post low and post high mixed blocks into quartiles for div model
    % with smaller integration window
    [divModel.postLow(rr,:), divModel.postHigh(rr,:), divModel.postLow_q1(rr,:),  ...
        divModel.postHigh_q1(rr, :)] = quartileAnalysis_SS(div);
    
    % detrend expert rat wait times 
    E.A = detrendwt_SS(E.A);

    % conditional wait time for previous 10 rewards
        % Use noisy wait times for inference model to show potential spread
        % Predicted wait times for inference model are drawn from noise
        % distribution and will not be exactly the same each run
    inf.wait_time = inf_wt_mdl; 
    [expert.prevSum_norm(rr,:), expert.prevSum_raw(rr,:)] = ...
        prev10Rews(E.A, sumtrials, BlkInf); 
    infModel.prevSum_norm(rr,:) = prev10Rews(inf, sumtrials, BlkInf);
    divModel.prevSum_norm(rr,:) = prev10Rews(div, sumtrials, BlkInf);

    % conditional wait time for previous reward
    [expert.prevRew(rr,:), ~] = waittime_20ul_by_previous_vol_SS(E.A, BlkInf);
    [infModel.prevRew(rr,:), ~] = waittime_20ul_by_previous_vol_SS(inf, BlkInf);
    [divModel.prevRew(rr,:), ~] = waittime_20ul_by_previous_vol_SS(div, BlkInf);

end

%% Averages

div_mtol = mean(divModel.mtol, 'omitnan');
div_mtol_sem = sem(divModel.mtol);

div_mtoh = mean(divModel.mtoh, 'omitnan');
div_mtoh_sem = sem(divModel.mtoh);

div_ltom = mean(divModel.ltom, 'omitnan');
div_ltom_sem = sem(divModel.ltom);

div_htom = mean(divModel.htom, 'omitnan');
div_htom_sem = sem(divModel.htom);

div_postLow = mean(divModel.postLow, 'omitnan');
div_postLow_sem = sem(divModel.postLow);

div_postHigh = mean(divModel.postHigh, 'omitnan');
div_postHigh_sem = sem(divModel.postHigh);

div_postLowq1 = mean(divModel.postLow_q1, 'omitnan');
div_postLowq1_sem = sem(divModel.postLow_q1);

div_postHighq1 = mean(divModel.postHigh_q1, 'omitnan');
div_postHighq1_sem = sem(divModel.postHigh_q1);

%% Conditional wait time stats

pExp = signrank(expert.prevSum_norm(:,1), expert.prevSum_norm(:,2)); 
pInf = signrank(infModel.prevSum_norm(:,1), infModel.prevSum_norm(:,2));
pDiv = signrank(divModel.prevSum_norm(:,1), divModel.prevSum_norm(:,2));

pExp_rew = signrank(expert.prevRew(:,1), expert.prevRew(:,2)); 
pInf_rew = signrank(infModel.prevRew(:,1), infModel.prevRew(:,2));
pDiv_rew = signrank(divModel.prevRew(:,1), divModel.prevRew(:,2));

%Effect sizes
sumEffect_exp = mean(expert.prevSum_raw(:,1), 'omitnan') - ...
    mean(expert.prevSum_raw(:,2), 'omitnan'); 

%% Fit sigmoids to simulated transition dynamics for 10 div model and inference model
disp('Fitting sigmoids to simulations')
disp('This will take ~ 5 minutes')

wndw = [find(xvec==-10):find(xvec==25)];

% rng(20)
% useRats = randsample(nrats, 50);

%rat trial data used:
useRats = [274; 93; 220; 161; 263; 269; 285; 242; 230; 206; 257; 133; 228;...
238; 181; 115; 13; 314; 85; 112; 91; 89; 293; 197; 138; 300; 57; 209; 96; ...
168; 68; 204; 63; 221; 195; 94; 270; 312; 109; 24; 55; 41; 332; 271; 179; ...
289; 174; 251; 297; 172];

n = length(useRats);

finalParams_ltom_inf = nan(n, 4);
finalParams_mtol_inf = nan(n, 4);
finalParams_htom_inf = nan(n, 4);
finalParams_mtoh_inf = nan(n, 4);
finalParams_ltom_inf_linear = nan(n, 4);
finalParams_mtol_inf_linear = nan(n, 4);
finalParams_htom_inf_linear = nan(n, 4);
finalParams_mtoh_inf_linear = nan(n, 4);

finalParams_ltom_div = nan(n, 4);
finalParams_mtol_div = nan(n, 4);
finalParams_htom_div = nan(n, 4);
finalParams_mtoh_div = nan(n, 4);
finalParams_ltom_div_linear = nan(n, 4);
finalParams_mtol_div_linear = nan(n, 4);
finalParams_htom_div_linear = nan(n, 4);
finalParams_mtoh_div_linear = nan(n, 4);

for rr = 1:length(useRats)
    try %skip if there are nans in the transition dynamics data due to no catch trial at that bin - cant fit sigmoid
        %inference model 
            % log2 rewards
        finalParams_ltom_inf(rr, :) = my_fit_sigmoid(xvec(wndw), ...
            infModel.ltom(useRats(rr), wndw), 50);
        finalParams_mtol_inf(rr, :) = my_fit_sigmoid(xvec(wndw), ...
            infModel.mtol(useRats(rr), wndw), 50);
        finalParams_htom_inf(rr, :) = my_fit_sigmoid(xvec(wndw), ...
            infModel.htom(useRats(rr), wndw), 50);
        finalParams_mtoh_inf(rr, :) = my_fit_sigmoid(xvec(wndw), ...
            infModel.mtoh(useRats(rr), wndw), 50);
            % linear rewards
        finalParams_ltom_inf_linear(rr, :) = my_fit_sigmoid(xvec(wndw), ...
            infModel.ltom_linear(useRats(rr), wndw), 50);
        finalParams_mtol_inf_linear(rr, :) = my_fit_sigmoid(xvec(wndw), ...
            infModel.mtol_linear(useRats(rr), wndw), 50);
        finalParams_htom_inf_linear(rr, :) = my_fit_sigmoid(xvec(wndw), ...
            infModel.htom_linear(useRats(rr), wndw), 50);
        finalParams_mtoh_inf_linear(rr, :) = my_fit_sigmoid(xvec(wndw), ...
            infModel.mtoh_linear(useRats(rr), wndw), 50);

        % divisive normalization model with smaller integration window
            % log 2 rewards
        finalParams_ltom_div(rr, :) = my_fit_sigmoid(xvec(wndw), ...
            divModel.ltom(useRats(rr), wndw), 50);
        finalParams_mtol_div(rr, :) = my_fit_sigmoid(xvec(wndw), ...
            divModel.mtol(useRats(rr), wndw), 50);
        finalParams_htom_div(rr, :) = my_fit_sigmoid(xvec(wndw), ...
            divModel.htom(useRats(rr), wndw), 50);
        finalParams_mtoh_div(rr, :) = my_fit_sigmoid(xvec(wndw), ...
            divModel.mtoh(useRats(rr), wndw), 50);
            % linear rewards
        finalParams_ltom_div_linear(rr, :) = my_fit_sigmoid(xvec(wndw), ...
            divModel.ltom_linear(useRats(rr), wndw), 50);
        finalParams_mtol_div_linear(rr, :) = my_fit_sigmoid(xvec(wndw), ...
            divModel.mtol_linear(useRats(rr), wndw), 50);
        finalParams_htom_div_linear(rr, :) = my_fit_sigmoid(xvec(wndw), ...
            divModel.htom_linear(useRats(rr), wndw), 50);
        finalParams_mtoh_div_linear(rr, :) = my_fit_sigmoid(xvec(wndw), ...
            divModel.mtoh_linear(useRats(rr), wndw), 50);
    catch
        continue
    end
end

%% medians and stats of sigmoid fits
i_mtol = median(finalParams_mtol_inf(:,3), 'omitnan');
i_ltom = median(finalParams_ltom_inf(:,3), 'omitnan');
i_mtoh = median(finalParams_mtoh_inf(:,3), 'omitnan'); 
i_htom = median(finalParams_htom_inf(:,3), 'omitnan');

i_mtol_linear = median(finalParams_mtol_inf_linear(:,3), 'omitnan');
i_ltom_linear = median(finalParams_ltom_inf_linear(:,3), 'omitnan');
i_mtoh_linear = median(finalParams_mtoh_inf_linear(:,3), 'omitnan'); 
i_htom_linear = median(finalParams_htom_inf_linear(:,3), 'omitnan');

pLow_inf = signrank(finalParams_mtol_inf(:,3), finalParams_ltom_inf(:,3));
pHigh_inf = signrank(finalParams_mtoh_inf(:,3), finalParams_htom_inf(:,3));
pLow_inf_linear = signrank(finalParams_mtol_inf_linear(:,3), finalParams_ltom_inf_linear(:,3));
pHigh_inf_linear = signrank(finalParams_mtoh_inf_linear(:,3), finalParams_htom_inf_linear(:,3));


d_mtol = median(finalParams_mtol_div(:,3), 'omitnan');
d_ltom = median(finalParams_ltom_div(:,3), 'omitnan');
d_mtoh = median(finalParams_mtoh_div(:,3), 'omitnan'); 
d_htom = median(finalParams_htom_div(:,3), 'omitnan');

d_mtol_linear = median(finalParams_mtol_div_linear(:,3), 'omitnan');
d_ltom_linear = median(finalParams_ltom_div_linear(:,3), 'omitnan');
d_mtoh_linear = median(finalParams_mtoh_div_linear(:,3), 'omitnan'); 
d_htom_linear = median(finalParams_htom_div_linear(:,3), 'omitnan');

pLow_div = signrank(finalParams_mtol_div(:,3), finalParams_ltom_div(:,3));
pHigh_div = signrank(finalParams_mtoh_div(:,3), finalParams_htom_div(:,3));
pLow_div_linear = signrank(finalParams_mtol_div_linear(:,3), finalParams_ltom_div_linear(:,3));
pHigh_div_linear = signrank(finalParams_mtoh_div_linear(:,3), finalParams_htom_div_linear(:,3));

%% Plot trainsition dynamics and quartile analysis for 10 trial divisive normalizaiton model

figure; hold on
tiledlayout(1, 4, 'TileSpacing', 'compact')

% Panel A - plot wait time dynamics for divisive normalization model with 10 trial window
nexttile
shadedErrorBar(xvec, div_mtol, div_mtol_sem, ...
    'lineprops', {'b', 'linewidth', 1});
shadedErrorBar(xvec, div_mtoh, div_mtoh_sem, ...
    'lineprops', {'r', 'linewidth', 1});
ylim([-1.5 1.5])
yl = ylim;
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('Wait time (z-scored)');
axis square
ax1 = gca;
ax1.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off
title('\rm Mixed into adaptation')

nexttile
shadedErrorBar(xvec, div_ltom, div_ltom_sem, ...
    'lineprops', {'b', 'linewidth', 1});
shadedErrorBar(xvec, div_htom, div_htom_sem, ...
    'lineprops', {'r', 'linewidth', 1});
ylim([-1.5 1.5])
yl = ylim;
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('Wait time (z-scored)');
axis square
ax2 = gca;
ax2.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off
title('\rm Adaptation into mixed')

% Plot quartile analysis
nexttile
shadedErrorBar(q, div_postLow, div_postLow_sem, 'lineprops', {'color', ...
    [0.1 0.1 0.6], 'linewidth', 1})
hold on
shadedErrorBar(q, div_postHigh, div_postHigh_sem, 'lineprops', {'color', ...
    [0.6 0.1 0.1], 'linewidth', 1})
set(gca, 'TickDir', 'out'); box off;
xlim([0 5]);
ylim([-0.5 0.5])
xticks([1:4])
xlabel('Quartile')
ylabel('Wait time (z-scored)')
ax3 = gca;
ax3.YRuler.TickLabelGapOffset = 1;
axis square

nexttile
shadedErrorBar(1:5, div_postLowq1, div_postLowq1_sem, 'lineprops', {'color', ...
    [0.1 0.1 0.6], 'linewidth', 1})
hold on
shadedErrorBar(1:5, div_postHighq1, div_postHighq1_sem, 'lineprops', {'color', ...
    [0.6 0.1 0.1], 'linewidth', 1})
set(gca, 'xTick', 1:5);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
ylim([10 45])
xlabel('Reward offer')
ylabel('Wait time (s)')
ax4 = gca;
ax4.YRuler.TickLabelGapOffset = 1;
axis square
set(gcf,'renderer','painter')
set(gcf, 'units', 'centimeters', 'position', [10 10 20 4])

%% Plot conditional wait time for 20ul panels

x = [0.5 1];
% avgDiff = mean(expert.prevSum_raw(:,1)) - mean(expert.prevSum_raw(:,2));

figure; hold on
tiledlayout(1, 6, 'TileSpacing', 'compact')

nexttile
hold on
arrayfun(@(k) plot(x, expert.prevRew(k,:) , 'color', [0.65 0.65 0.65 0.4]), 1:nrats)
plot(x, mean(expert.prevRew, 'omitnan'), 'k', 'LineWidth', 2)
xlim([0.4 1.1])
ylim([-1 2.2])
xticks([0.5 1])
xticklabels({'<20' '>20'})
xlabel('Previous reward')
ylabel('Wait time (z-scored)')
set(gca, 'TickDir', 'out'); box off;
ax1 = gca;
ax1.YRuler.TickLabelGapOffset = 1;
% axis square
title('\rm Expert')

nexttile;
hold on
arrayfun(@(k) plot(x, infModel.prevRew(k,:) , 'color', [0.65 0.65 0.65 0.4]), 1:nrats)
plot(x, mean(infModel.prevRew, 'omitnan'), 'k', 'LineWidth', 2)
xlim([0.4 1.1])
ylim([-1 2.2])
xticks([0.5 1])
xticklabels({'<20' '>20'})
xlabel('Previous reward')
ylabel('Wait time (z-scored)')
set(gca, 'TickDir', 'out'); box off;
ax2 = gca;
ax2.YRuler.TickLabelGapOffset = 1;
% axis square
title('\rm Inference model')

nexttile;
hold on
arrayfun(@(k) plot(x, divModel.prevRew(k,:) , 'color', [0.65 0.65 0.65 0.4]), 1:nrats)
plot(x, mean(divModel.prevRew, 'omitnan'), 'k', 'LineWidth', 2)
xlim([0.4 1.1])
ylim([-1 2.2])
xticks([0.5 1])
xticklabels({'<20' '>20'})
xlabel('Previous reward')
ylabel('Wait time (z-scored)')
set(gca, 'TickDir', 'out'); box off;
ax3 = gca;
ax3.YRuler.TickLabelGapOffset = 1;
% axis square
title('\rm Div norm model')


nexttile;
hold on
arrayfun(@(k) plot(x, expert.prevSum_norm(k,:) , 'color', [0.65 0.65 0.65 0.4]), 1:nrats)
plot(x, mean(expert.prevSum_norm, 'omitnan'), 'k', 'LineWidth', 2)
xlim([0.4 1.1])
ylim([-1 1])
xticks([0.5 1])
xticklabels({'50%' '50%'})
xlabel('Previous sum')
ylabel('Wait time (z-scored)')
set(gca, 'TickDir', 'out'); box off;
ax4 = gca;
ax4.YRuler.TickLabelGapOffset = 1;
% axis square
title('\rm Expert')

nexttile;
hold on
arrayfun(@(k) plot(x, infModel.prevSum_norm(k,:) , 'color', [0.65 0.65 0.65 0.4]), 1:nrats)
plot(x, mean(infModel.prevSum_norm, 'omitnan'), 'k', 'LineWidth', 2)
xlim([0.4 1.1])
ylim([-1 1])
xticks([0.5 1])
xticklabels({'50%' '50%'})
xlabel('Previous sum')
ylabel('Wait time (z-scored)')
set(gca, 'TickDir', 'out'); box off;
ax5 = gca;
ax5.YRuler.TickLabelGapOffset = 1;
% axis square
title('\rm Inference model')

nexttile;
hold on
arrayfun(@(k) plot(x, divModel.prevSum_norm(k,:) , 'color', [0.65 0.65 0.65 0.4]), 1:nrats)
plot(x, mean(divModel.prevSum_norm, 'omitnan'), 'k', 'LineWidth', 2)
xlim([0.4 1.1])
ylim([-1 1])
xticks([0.5 1])
xticklabels({'50%' '50%'})
xlabel('Previous sum')
ylabel('Wait time (z-scored)')
set(gca, 'TickDir', 'out'); box off;
ax6 = gca;
ax6.YRuler.TickLabelGapOffset = 1;
% axis square
title('\rm Div norm model')


set(gcf,'renderer','painter')
set(gcf, 'units', 'centimeters', 'position', [10 10 20 4])

%% 
figure;
tiledlayout(2, 4, 'tilespacing', 'compact')

%HISTOGRAMS
nexttile
histogram(finalParams_mtol_inf(:,3), 10,'facecolor', 'b')
ylim([0 20])
xlim([0 10])
hold on
xline(i_mtol, '--b')
histogram(finalParams_ltom_inf(:,3), 10,'facecolor', 'k')
xline(i_ltom, '--k')
set(gca, 'TickDir', 'out'); box off;
xlabel('X offset')
ylabel('# simulations')

nexttile
histogram(finalParams_mtoh_inf(:,3), 10, 'facecolor', 'r')
ylim([0 20])
xlim([0 10])
hold on
xline(i_mtoh, '--r')
histogram(finalParams_htom_inf(:,3), 10, 'facecolor', 'k')
xline(i_htom, '--k')
set(gca, 'TickDir', 'out'); box off;
xlabel('X offset')
ylabel('# simulations')

nexttile
histogram(finalParams_mtol_inf_linear(:,3), 10,'facecolor', 'b')
ylim([0 20])
xlim([0 10])
hold on
xline(i_mtol, '--b')
histogram(finalParams_ltom_inf_linear(:,3), 10,'facecolor', 'k')
xline(i_ltom, '--k')
set(gca, 'TickDir', 'out'); box off;
xlabel('X offset')
ylabel('# simulations')

nexttile
histogram(finalParams_mtoh_inf_linear(:,3), 10, 'facecolor', 'r')
ylim([0 20])
xlim([0 10])
hold on
xline(i_mtoh, '--r')
histogram(finalParams_htom_inf_linear(:,3), 10, 'facecolor', 'k')
xline(i_htom, '--k')
set(gca, 'TickDir', 'out'); box off;
xlabel('X offset')
ylabel('# simulations')

nexttile
histogram(finalParams_mtol_div(:,3), 10,'facecolor', 'b')
ylim([0 20])
xlim([0 10])
hold on
xline(d_mtol, '--b')
histogram(finalParams_ltom_div(:,3), 10,'facecolor', 'k')
xline(d_ltom, '--k')
set(gca, 'TickDir', 'out'); box off;
xlabel('X offset')
ylabel('# simulations')

nexttile
histogram(finalParams_mtoh_div(:,3), 10, 'facecolor', 'r')
ylim([0 20])
xlim([0 10])
hold on
xline(d_mtoh, '--r')
histogram(finalParams_htom_div(:,3), 10, 'facecolor', 'k')
xline(d_htom, '--k')
set(gca, 'TickDir', 'out'); box off;
xlabel('X offset')
ylabel('# simulations')

nexttile
histogram(finalParams_mtol_div_linear(:,3), 10,'facecolor', 'b')
ylim([0 20])
xlim([0 10])
hold on
xline(d_mtol, '--b')
histogram(finalParams_ltom_div_linear(:,3), 10,'facecolor', 'k')
xline(d_ltom, '--k')
set(gca, 'TickDir', 'out'); box off;
xlabel('X offset')
ylabel('# simulations')

nexttile
histogram(finalParams_mtoh_div_linear(:,3), 10, 'facecolor', 'r')
ylim([0 20])
xlim([0 10])
hold on
xline(d_mtoh, '--r')
histogram(finalParams_htom_div_linear(:,3), 10, 'facecolor', 'k')
xline(d_htom, '--k')
set(gca, 'TickDir', 'out'); box off;
xlabel('X offset')
ylabel('# simulations')

set(gcf,'renderer','painter')
set(gcf, 'Color', [1 1 1]);
set(gcf, 'units', 'centimeters', 'position', [10 10 20 8])
