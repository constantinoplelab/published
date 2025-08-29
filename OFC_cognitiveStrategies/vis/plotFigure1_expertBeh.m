function plotFigure1_expertBeh(processedBehaviorPath)
%Plot figure 1 - expert rat behavior and inference model predictions
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

e_mtol = mean(expert.mtol, 'omitnan');
e_mtol_sem = sem(expert.mtol);

e_mtoh = mean(expert.mtoh, 'omitnan');
e_mtoh_sem = sem(expert.mtoh);

e_ltom = mean(expert.ltom, 'omitnan');
e_ltom_sem = sem(expert.ltom);

e_htom = mean(expert.htom, 'omitnan');
e_htom_sem = sem(expert.htom);

inf_mtol = mean(infModel.mtol, 'omitnan');
inf_mtol_sem = sem(infModel.mtol);

inf_mtoh = mean(infModel.mtoh, 'omitnan');
inf_mtoh_sem = sem(infModel.mtoh);

inf_ltom = mean(infModel.ltom, 'omitnan');
inf_ltom_sem = sem(infModel.ltom);

inf_htom = mean(infModel.htom, 'omitnan');
inf_htom_sem = sem(infModel.htom);

e_postLow = mean(expert.postLow, 'omitnan');
e_postLow_sem = sem(expert.postLow);

e_postHigh = mean(expert.postHigh, 'omitnan');
e_postHigh_sem = sem(expert.postHigh);

e_postLowq1 = mean(expert.postLow_q1, 'omitnan'); %from incongruent
e_postLowq1_sem = sem(expert.postLow_q1);

e_postHighq1 = mean(expert.postHigh_q1, 'omitnan');
e_postHighq1_sem = sem(expert.postHigh_q1);

e_postLowq1_con = mean(expert.postLow_q1Con, 'omitnan'); %to first incongruent
e_postLowq1_con_sem = sem(expert.postLow_q1Con);

e_postHighq1_con = mean(expert.postHigh_q1Con, 'omitnan');
e_postHighq1_con_sem = sem(expert.postHigh_q1Con);


inf_postLow = mean(infModel.postLow, 'omitnan');
inf_postLow_sem = sem(infModel.postLow);

inf_postHigh = mean(infModel.postHigh, 'omitnan');
inf_postHigh_sem = sem(infModel.postHigh);

inf_postLowq1 = mean(infModel.postLow_q1, 'omitnan');
inf_postLowq1_sem = sem(infModel.postLow_q1);

inf_postHighq1 = mean(infModel.postHigh_q1, 'omitnan');
inf_postHighq1_sem = sem(infModel.postHigh_q1);

inf_postLowq1_con = mean(infModel.postLow_q1Con, 'omitnan');
inf_postLowq1_con_sem = sem(infModel.postLow_q1Con);

inf_postHighq1_con = mean(infModel.postHigh_q1Con, 'omitnan');
inf_postHighq1_con_sem = sem(infModel.postHigh_q1Con);


loAvg = mean(expert.lo(:, 1:3), 2, 'omitnan');
loConAvg = mean(expert.postLow_q1Con(:,1:3), 2, 'omitnan');
hiAvg = mean(expert.hi(:, 3:5), 2, 'omitnan');
hiConAvg = mean(expert.postHigh_q1Con(:,3:5), 2, 'omitnan');

%% stats

p_expwt = signrank(expert.lo(:,3), expert.hi(:,3));
p_expL = anova1(expert.postLow(:, 2:end), [], 'off');
p_expH = anova1(expert.postHigh(:, 2:end), [], 'off');
p_expQ1 = (arrayfun(@(x) signrank(expert.postLow_q1(:,x), ...
    expert.postHigh_q1(:,x)), 1:5)).*5;
p_expQ1Con = signrank(expert.postLow_q1Con(:,3), expert.postHigh_q1Con(:,3));
p_expConLow = signrank(loAvg, loConAvg);
p_expConHigh = signrank(hiAvg, hiConAvg);

%% fit sigmoids to individual rat transition dynamics
disp('Fitting sigmoid functions to transition dynamics')
disp('This will take ~5-7 minutes')

%window to fit sigmoid -- produces the most reasonable fits
wndw = [find(xvec==-10):find(xvec==25)];

%aligned to true block transition 
nrats = length(expert.ltom);
finalParams_ltom_exp = nan(nrats, 4);
finalParams_mtol_exp = nan(nrats, 4);
finalParams_htom_exp = nan(nrats, 4);
finalParams_mtoh_exp = nan(nrats, 4);

for rr = 1:nrats
    try
        finalParams_ltom_exp(rr, :) = my_fit_sigmoid(xvec(wndw), expert.ltom(rr, wndw), 50);
        finalParams_mtol_exp(rr, :) = my_fit_sigmoid(xvec(wndw), expert.mtol(rr, wndw), 50);
        finalParams_htom_exp(rr, :) = my_fit_sigmoid(xvec(wndw), expert.htom(rr, wndw), 50);
        finalParams_mtoh_exp(rr, :) = my_fit_sigmoid(xvec(wndw), expert.mtoh(rr, wndw), 50);
    catch
        continue
    end
end

%aligned to incongruent 
finalParams_ltomIncong_exp = nan(nrats, 4);
finalParams_htomIncong_exp = nan(nrats, 4);
for rr = 1:nrats
    try
        finalParams_ltomIncong_exp(rr, :) = my_fit_sigmoid(xvec(wndw), ...
            expert.ltom_incong(rr, wndw), 50);  
        finalParams_htomIncong_exp(rr, :) = my_fit_sigmoid(xvec(wndw), ...
            expert.htom_incong(rr, wndw), 50);
    catch
        continue
    end
end

%% plot example block durations
figure; hold on
tiledlayout(2, 1, 'tilespacing', 'compact')

for ii = 1:2
    nexttile
    histogram(expert.blockDuration{ii}, 20, 'FaceColor', 'k'); 
    xlim([40 100]);
    set(gca, 'xTick', [40 60 80 100]);
    xlabel('# of trials');
    ylabel('N (Blocks)');
    set(gca, 'TickDir', 'out'); box off;
end
set(gcf, 'units', 'centimeters', 'position', [10 10 4 4])
set(gcf,'renderer','painter')

%% plot expert behavior and inference model predictions

figure; hold on
tiledlayout(2, 6, 'TileSpacing', 'compact')

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
ylim([9 13.2])
xlabel('Reward offer')
ylabel('Mean wait time (s)')
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
shadedErrorBar(xvec, e_ltom, e_ltom_sem, ...
    'lineprops', {'-b', 'linewidth', 1}); hold on
shadedErrorBar(xvec, e_htom, e_htom_sem, ...
    'lineprops', {'-r', 'linewidth', 1});
ylim([-0.2 0.3])
yl = ylim;
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
yticks([-0.2, 0, 0.2])
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
ylim([-0.12 0.1])
xticks([1:4])
yticks([-0.1, 0, 0.1])
xlabel('Quartile')
ylabel({'Mean z-scored', 'wait time'})
text(1, 0.14, 'Post-low', 'color', [0.1 0.1 0.6], 'FontSize', 8)
text(1, 0.11, 'Post-high', 'color', [0.6 0.1 0.1], 'FontSize', 8)
ax4 = gca;
ax4.YRuler.TickLabelGapOffset = 1;
% axis square

nexttile
shadedErrorBar(x, e_postLowq1_con, e_postLowq1_con_sem, 'lineprops', {'color', ...
    [0.1 0.1 0.6], 'linewidth', 1})
hold on
shadedErrorBar(x, e_postHighq1_con, e_postHighq1_con_sem, 'lineprops', {'color', ...
    [0.6 0.1 0.1], 'linewidth', 1})
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
ylim([9 13.2])
xlabel('Reward offer')
ylabel('Wait time (s)')
ax5 = gca;
ax5.YRuler.TickLabelGapOffset = 1;
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
ylim([9 13.2])
xlabel('Reward offer')
ylabel('Wait time (s)')
ax6 = gca;
ax6.YRuler.TickLabelGapOffset = 1;
% axis square

nexttile %leave blank for model schematics


%wait time dynamics -- inference model
nexttile
shadedErrorBar(xvec, inf_mtol, inf_mtol_sem, 'lineprops', ...
    {'-b', 'linewidth', 1})
hold on
shadedErrorBar(xvec, inf_mtoh, inf_mtoh_sem, ...
    'lineprops', {'-r', 'linewidth', 1});
ylim([-2 2])
yl = ylim;
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('Wait time (z-scored)');
% axis square
ax7 = gca;
ax7.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off

nexttile
shadedErrorBar(xvec, inf_ltom, inf_ltom_sem, ...
    'lineprops', {'-b', 'linewidth', 1}); hold on
shadedErrorBar(xvec, inf_htom, inf_htom_sem, ...
    'lineprops', {'-r', 'linewidth', 1});
ylim([-2 2])
yl = ylim;
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
% ylabel('\Delta z-scored wait time');
% axis square
ax8 = gca;
ax8.YRuler.TickLabelGapOffset = 1;
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
ylim([-0.3 0.3])
xticks([1:4])
xlabel('Quartile')
ylabel({'Mean z-scored', 'wait time'})
% text(1, 0.5, 'Post-low', 'color', [0.1 0.1 0.6], 'FontSize', 8)
% text(1, 0.4, 'Post-high', 'color', [0.6 0.1 0.1], 'FontSize', 8)
ax9 = gca;
ax9.YRuler.TickLabelGapOffset = 1;
% axis square

nexttile
shadedErrorBar(x, inf_postLowq1_con, inf_postLowq1_con_sem, 'lineprops', {'color', ...
    [0.1 0.1 0.6], 'linewidth', 1})
hold on
shadedErrorBar(x, inf_postHighq1_con, inf_postHighq1_con_sem, 'lineprops', {'color', ...
    [0.6 0.1 0.1], 'linewidth', 1})
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
ylim([9 13.5])
xlabel('Reward offer')
ylabel('Wait time (s)')
ax10 = gca;
ax10.YRuler.TickLabelGapOffset = 1;

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
ylim([9 13.5])
xlabel('Reward offer')
ylabel('Wait time (s)')
ax11 = gca;
ax11.YRuler.TickLabelGapOffset = 1;
% axis square

set(gcf,'renderer','painter')
set(gcf, 'Color', [1 1 1]);
set(gcf, 'units', 'centimeters', 'position', fsize)


%% plot histograms of sigmoid fits to each rat -- expert
e_mtol = median(finalParams_mtol_exp(:,3), 'omitnan');
e_ltom = median(finalParams_ltom_exp(:,3), 'omitnan');
e_mtoh = median(finalParams_mtoh_exp(:,3), 'omitnan'); 
e_htom = median(finalParams_htom_exp(:,3), 'omitnan');
e_ltomIncong = median(finalParams_ltomIncong_exp(:,3), 'omitnan');
e_htomIncong = median(finalParams_htomIncong_exp(:,3), 'omitnan');

p_low = signrank(finalParams_mtol_exp(:,3), finalParams_ltom_exp(:,3));
p_high = signrank(finalParams_mtoh_exp(:,3), finalParams_htom_exp(:,3));

p_incong = signrank(finalParams_ltomIncong_exp(:,3), finalParams_htomIncong_exp(:,3));
p_trueChange = signrank(finalParams_ltom_exp(:,3), finalParams_htom_exp(:,3));

%example sigmoid fits
ex1_ltom = finalParams_ltom_exp(302,4) + (finalParams_ltom_exp(302,1) ...
    -finalParams_ltom_exp(302,4)) ./ (1 + exp(-finalParams_ltom_exp(302,2)*...
    (xvec-finalParams_ltom_exp(302,3))));
ex1_mtol = finalParams_mtol_exp(302,4) + (finalParams_mtol_exp(302,1) ...
    -finalParams_mtol_exp(302,4)) ./ (1 + exp(-finalParams_mtol_exp(302,2)*...
    (xvec-finalParams_mtol_exp(302,3))));
ex1_htom = finalParams_htom_exp(302,4) + (finalParams_htom_exp(302,1) ...
    -finalParams_htom_exp(302,4)) ./ (1 + exp(-finalParams_htom_exp(302,2)*...
    (xvec-finalParams_htom_exp(302,3))));
ex1_mtoh = finalParams_mtoh_exp(302,4) + (finalParams_mtoh_exp(302,1) ...
    -finalParams_mtoh_exp(302,4)) ./ (1 + exp(-finalParams_mtoh_exp(302,2)*...
    (xvec-finalParams_mtoh_exp(302,3))));

%Sem order: ltom, htom, mtol, mtoh, ltom_incong, htom_incong

figure; hold on
tiledlayout(1, 6, 'TileSpacing', 'compact');

%EXAMPLE RAT 1
nexttile
shadedErrorBar(xvec, expert.mtol(302,:), expert.dynamicsSems{302}(3,:), ...
    'lineprops', {'-b', 'linewidth', 1})
hold on
shadedErrorBar(xvec, expert.mtoh(302,:), expert.dynamicsSems{302}(4,:), ...
    'lineprops', {'-r', 'linewidth', 1});
plot(xvec, ex1_mtol, '--k', 'linewidth', 1)
plot(xvec, ex1_mtoh, '--k', 'linewidth', 1)
ylim([-0.5 0.7])
yl = ylim;
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Trial from block switch');
ylabel('Wait time (z-scored)');
% axis square
ax1 = gca;
ax1.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off
title('Example rat')

nexttile
shadedErrorBar(xvec, expert.ltom(302,:), expert.dynamicsSems{302}(1,:), ...
    'lineprops', {'-b', 'linewidth', 1}); hold on
shadedErrorBar(xvec, expert.htom(302,:), expert.dynamicsSems{302}(2,:), ...
    'lineprops', {'-r', 'linewidth', 1});
plot(xvec, ex1_ltom, '--k', 'linewidth', 1)
plot(xvec, ex1_htom, '--k', 'linewidth', 1)
ylim([-0.5 0.7])
yl = ylim;
xlim([-25 40]);
line([0 0], [yl(1) yl(2)], 'Color', [0 0 0], 'LineStyle', '--');
% axis square
ax2 = gca;
ax2.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off


%HISTOGRAMS
nexttile
histogram(finalParams_mtol_exp(:,3), 'facecolor', 'b')
yl = ylim;
hold on
xline(e_mtol, '--b')
histogram(finalParams_ltom_exp(:,3), 'facecolor', 'k')
xline(e_ltom, '--k')
set(gca, 'TickDir', 'out'); box off;
xlabel('X offset')
ylabel('# rats')
title(['\rm p = ' string(p_low)])


nexttile
histogram(finalParams_mtoh_exp(:,3), 'facecolor', 'r')
yl = ylim;
hold on
xline(e_mtoh, '--r')
histogram(finalParams_htom_exp(:,3), 'facecolor', 'k')
xline(e_htom, '--k')
set(gca, 'TickDir', 'out'); box off;
xlabel('X offset')
ylabel('# rats')
title(['\rm p = ' string(p_high)])

nexttile
hold on
histogram(finalParams_ltomIncong_exp(:,3), 'facecolor', 'k')
xline(e_ltomIncong, '--k')
set(gca, 'TickDir', 'out'); box off;
xlabel('X offset')
ylabel('# rats')
title(['\rm' string(e_ltomIncong)])

nexttile
hold on
histogram(finalParams_htomIncong_exp(:,3), 'facecolor', 'k')
xline(e_htomIncong, '--k')
set(gca, 'TickDir', 'out'); box off;
xlabel('X offset')
ylabel('# rats')
title(['\rm' string(e_htomIncong)])

set(gcf,'renderer','painter')
set(gcf, 'units', 'centimeters', 'position', [10 10 20 4)


