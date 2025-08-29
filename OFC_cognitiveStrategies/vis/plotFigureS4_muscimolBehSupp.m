function plotFigureS4_muscimolBehSupp(processedMuscimolPath)
% Plot Figure S4. YOU MUST RUN processMuscimolData BEFORE RUNNING THIS
% FUNCTION

% INPUTS:
%   processedMuscimolPath = local path to output saved from running
%       processMuscimolData

%% 
x = 1:5;

%general functions
sem = @(x) std(x,'omitnan')./sqrt(sum(any(~isnan(x), 2)));
setDefaultFigProps

%% load processed muscimol behavior data
load([processedMuscimolPath, filesep, 'control.mat'])
load([processedMuscimolPath, filesep, 'muscimol.mat'])

nrats = length(control.mix);

%% averages and p values

%normalized wait times
normMix_C = cell2mat(arrayfun(@(x) control.mix(x,:)./control.mix(x, 3), ...
    1:nrats, 'uniformoutput', false)');
normHi_C = cell2mat(arrayfun(@(x) control.hi(x,:)./control.mix(x, 3), ...
    1:nrats, 'uniformoutput', false)');
normLo_C = cell2mat(arrayfun(@(x) control.lo(x,:)./control.mix(x, 3), ...
    1:nrats, 'uniformoutput', false)');

normMix_M = cell2mat(arrayfun(@(x) muscimol.mix(x,:)./control.mix(x, 3), ...
    1:nrats, 'uniformoutput', false)');
normHi_M = cell2mat(arrayfun(@(x) muscimol.hi(x,:)./control.mix(x, 3), ...
    1:nrats, 'uniformoutput', false)');
normLo_M = cell2mat(arrayfun(@(x) muscimol.lo(x,:)./control.mix(x, 3), ...
    1:nrats, 'uniformoutput', false)');


%delta wait time ratios
wt1 = control.hi(:,3)./control.lo(:,3);
wt2 = muscimol.hi(:,3)./muscimol.lo(:,3);

deltawt = wt2 - wt1;

bins = -0.5:0.025:0.5;

pratio = signrank(wt2, wt1);


%slopes
pslopes_mix = signrank(control.slopes.mix(:,2), muscimol.slopes.mix(:,2));
pslopes_low = signrank(control.slopes.lo(:,2), muscimol.slopes.lo(:,2)); 
pslopes_high = signrank(control.slopes.hi(:,2), muscimol.slopes.hi(:,2));


%offets
poffsets_mix = signrank(control.slopes.mix(:,1), muscimol.slopes.mix(:,1));
poffsets_low = signrank(control.slopes.lo(:,1), muscimol.slopes.lo(:,1)); 
poffsets_high = signrank(control.slopes.hi(:,1), muscimol.slopes.hi(:,1));


%post-adapt q1 after first incongruent
postLowq1_C = mean(control.postLow_q1, 'omitnan');
postLowq1_sem_C = sem(control.postLow_q1);

postHighq1_C = mean(control.postHigh_q1, 'omitnan');
postHighq1_sem_C = sem(control.postHigh_q1);

postLowq1_M = mean(muscimol.postLow_q1, 'omitnan');
postLowq1_sem_M = sem(muscimol.postLow_q1);

postHighq1_M = mean(muscimol.postHigh_q1, 'omitnan');
postHighq1_sem_M = sem(muscimol.postHigh_q1);

p_conQ1 = (arrayfun(@(x) signrank(control.postLow_q1(:,x), ...
    control.postHigh_q1(:,x)), 1:5));
p_muscQ1 = (arrayfun(@(x) signrank(muscimol.postLow_q1(:,x), ...
    muscimol.postHigh_q1(:,x)), 1:5));


%ITIs
p_itiAll = signrank(control.all_iti.raw, muscimol.all_iti.raw);

iti_Cm = mean(control.mix_iti.z);
iti_Cm_sem = std(control.mix_iti.z);

iti_Ch = mean(control.hi_iti.z);
iti_Ch_sem = std(control.hi_iti.z);

iti_Cl = mean(control.lo_iti.z);
iti_Cl_sem = std(control.lo_iti.z);


iti_Mm = mean(muscimol.mix_iti.z);
iti_Mm_sem = std(muscimol.mix_iti.z);

iti_Mh = mean(muscimol.hi_iti.z);
iti_Mh_sem = std(muscimol.hi_iti.z);

iti_Ml = mean(muscimol.lo_iti.z);
iti_Ml_sem = std(muscimol.lo_iti.z);

% delta iti ratios
iti1 = control.hi_iti.raw ./ control.lo_iti.raw;
iti2 = muscimol.hi_iti.raw ./ muscimol.lo_iti.raw;

deltaiti = iti2 - iti1;

pratio_iti = signrank(iti2, iti1);

%% plot wait time curves and ratio

figure; hold on
tiledlayout(1, 3, 'TileSpacing', 'compact')

%wait time curves
nexttile
shadedErrorBar(x, mean(normHi_C), sem(normHi_C), 'lineprops', {'-r', ...
    'linewidth', 1.5})
shadedErrorBar(x, mean(normLo_C), sem(normLo_C), 'lineprops', {'-b', ...
    'linewidth', 1.5})
shadedErrorBar(x, mean(normMix_C), sem(normMix_C), 'lineprops', {'-k', ...
    'linewidth', 1.5})
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
ylim([0.7 1.5])
xlabel('Reward offer')
ylabel('Mean wait time (s)')
ax1 = gca;
ax1.YRuler.TickLabelGapOffset = 1;
% axis square
title('\rm Control')

nexttile
shadedErrorBar(x, mean(normHi_M), sem(normHi_M), 'lineprops', {'-r', ...
    'linewidth', 1.5})
shadedErrorBar(x, mean(normLo_M), sem(normLo_M), 'lineprops', {'-b', ...
    'linewidth', 1.5})
shadedErrorBar(x, mean(normMix_M), sem(normMix_M), 'lineprops', {'-k', ...
    'linewidth', 1.5})
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
ylim([0.7 1.5])
xlabel('Reward offer')
ylabel('Mean wait time (s)')
ax2 = gca;
ax2.YRuler.TickLabelGapOffset = 1;
% axis square
title('\rm Muscimol')

% delta wt histogram
nexttile
h = histogram(deltawt, bins);
h.FaceColor = [0.65 0.65 0.65];
h.EdgeColor = [0.65 0.65 0.65];
% text(-0.4, 2, strcat('p = ', num2str(pratio)))
xlim([-0.5 0.5])
ylim([0 3])
xline(0, 'k--', 'LineWidth', 1)
ylabel('N rats')    
xlabel({'\Delta wait time ratio'; '(muscimol - control)'})
% axis square
set(gca, 'TickDir', 'out'); box off;
ax3 = gca;
ax3.YRuler.TickLabelGapOffset = 1;
box off  

set(gcf,'renderer','painter')
set(gcf, 'units', 'centimeters', 'position', [10 10 16 4])

%% plot slopes and offsets

figure; hold on
tiledlayout(1, 2, 'tilespacing', 'compact')

%block slopes
jitx = randn(nrats, 6)*.03 ;
jitx = jitx + [0.25 0.75 1.75 2.25 3.25 3.75];

nexttile
scatter(jitx(:,1:2), [control.slopes.lo(:,2) ...
    muscimol.slopes.lo(:,2)], 15, 'b', 'filled')
hold on
arrayfun(@(k) plot(jitx(k,1:2), [control.slopes.lo(k,2) ...
    muscimol.slopes.lo(k,2)], 'b'), 1:nrats)
scatter(jitx(:,3:4), [control.slopes.mix(:,2) ...
    muscimol.slopes.mix(:,2)], 15, 'k', 'filled')
arrayfun(@(k) plot(jitx(k,3:4), [control.slopes.mix(k,2) ...
    muscimol.slopes.mix(k,2)], 'k'), 1:nrats)
scatter(jitx(:,5:6), [control.slopes.hi(:,2) ...
    muscimol.slopes.hi(:,2)], 15, 'r', 'filled')
arrayfun(@(k) plot(jitx(k,5:6), [control.slopes.hi(k,2) ...
    muscimol.slopes.hi(k,2)], 'r'), 1:nrats)
% text(1.5,1.3,strcat('p = ',num2str(pslopes_mix)))
xticks([0.25 0.75 1.75 2.25 3.25 3.75])
xlim([0 4])
ylim([-2.5 4.5])
ylabel('Slope')
xticklabels({'C' 'M' 'C' 'M' 'C' 'M'})
set(gca, 'TickDir', 'out'); box off;
ax4 = gca;
ax4.YRuler.TickLabelGapOffset = 1;
% axis square


%block offsets
nexttile
scatter(jitx(:,1:2), [control.slopes.lo(:,1) ...
    muscimol.slopes.lo(:,1)], 15, 'b', 'filled')
hold on
arrayfun(@(k) plot(jitx(k,1:2), [control.slopes.lo(k,1) ...
    muscimol.slopes.lo(k,1)], 'b'), 1:nrats)
scatter(jitx(:,3:4), [control.slopes.mix(:,1) ...
    muscimol.slopes.mix(:,1)], 15, 'k', 'filled')
arrayfun(@(k) plot(jitx(k,3:4), [control.slopes.mix(k,1) ...
    muscimol.slopes.mix(k,1)], 'k'), 1:nrats)
scatter(jitx(:,5:6), [control.slopes.hi(:,1) ...
    muscimol.slopes.hi(:,1)], 15, 'r', 'filled')
arrayfun(@(k) plot(jitx(k,5:6), [control.slopes.hi(k,1) ...
    muscimol.slopes.hi(k,1)], 'r'), 1:nrats)
% text(1.5,1.3,strcat('p = ',num2str(pslopes_mix)))
xlim([0 4])
ylim([-5 30])
ylabel('Offset')
xticks([0.25 0.75 1.75 2.25 3.25 3.75])
xticklabels({'C' 'M' 'C' 'M' 'C' 'M'})
set(gca, 'TickDir', 'out'); box off;
ax5 = gca;
ax5.YRuler.TickLabelGapOffset = 1;
% axis square

set(gcf,'renderer','painter')
set(gcf, 'units', 'centimeters', 'position', [10 10 16 4])


%% plot previous trial regressions and post-low/post-high 

figure; hold on
tiledlayout(1, 4, 'Tilespacing', 'Compact')

%regression plots 
nexttile
plot(0:6, control.regress(:, 1:7),...
    linewidth=0.5, color=[0.8 0.8 0.8])
shadedErrorBar(0:6, mean(control.regress(:,1:7)), sem(control.regress(:,1:7)),...
    lineprops={'color', 'k', 'linewidth', 2});
xticks(0:6)
xlim([-0.5, 6.5])
ylim([-0.1 0.6])
ax5 = gca;
ax5.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off;
% axis square
xlabel('Trials back')
ylabel('Regression coefficient')
title('\rm Control')

nexttile
plot(0:6, muscimol.regress(:, 1:7),...
    linewidth=0.5, color=[0.8 0.8 0.8])
shadedErrorBar(0:6, mean(muscimol.regress(:,1:7)), sem(muscimol.regress(:,1:7)),...
    lineprops={'color', 'k', 'linewidth', 2});
xticks(0:6)
xlim([-0.5, 6.5])
ylim([-0.1 0.6])
ax6 = gca;
ax6.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off;
% axis square
xlabel('Trials back')
ylabel('Regression coefficient')
title('\rm Muscimol')

%Post-low, post-high q1
nexttile
shadedErrorBar(x, postLowq1_C, postLowq1_sem_C, 'lineprops', {'color', ...
    [0.1 0.1 0.6], 'linewidth', 1})
hold on
shadedErrorBar(x, postHighq1_C, postHighq1_sem_C, 'lineprops', {'color', ...
    [0.6 0.1 0.1], 'linewidth', 1})
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
ylim([5 18])
xlabel('Reward offer')
ylabel('Wait time (s)')
ax5 = gca;
ax5.YRuler.TickLabelGapOffset = 1;
title('\rm Control')
% axis square

nexttile
shadedErrorBar(x, postLowq1_M, postLowq1_sem_M, 'lineprops', {'color', ...
    [0.1 0.1 0.6], 'linewidth', 1})
hold on
shadedErrorBar(x, postHighq1_M, postHighq1_sem_M, 'lineprops', {'color', ...
    [0.6 0.1 0.1], 'linewidth', 1})
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
ylim([5 18])
xlabel('Reward offer')
ylabel('Wait time (s)')
ax8 = gca;
ax8.YRuler.TickLabelGapOffset = 1;
title('\rm Muscimol')
% axis square

set(gcf,'renderer','painter')
set(gcf, 'units', 'centimeters', 'position', [10 10 16 4])


%% Plot ITIs

figure; hold on
tiledlayout(1, 3, 'Tilespacing', 'compact')

%ITIs
nexttile
errorbar(0.25, iti_Cl, iti_Cl_sem, 'bo', 'MarkerFaceColor', 'b', ...
    'LineWidth', 1, 'CapSize', 0)
hold on
errorbar(0.5, iti_Cm, iti_Cm_sem, 'ko', 'MarkerFaceColor', 'k', ...
    'LineWidth', 1, 'CapSize', 0)
errorbar(0.75, iti_Ch, iti_Ch_sem, 'ro', 'MarkerFaceColor', 'r', ...
    'LineWidth', 1, 'CapSize', 0)
xlim([0 1])
ylim([-0.2, 0.2])
xticks(0.25:0.25:0.75)
xticklabels({'Low', 'Mix', 'High'})
xlabel('Reward block')
ylabel('Mean trial initiation time (s)')
set(gca, 'TickDir', 'out'); box off;
% axis square
ax9 = gca;
ax9.YRuler.TickLabelGapOffset = 1;
title('\rm Control')

nexttile
errorbar(0.25, iti_Ml, iti_Ml_sem, 'bo', 'MarkerFaceColor', 'b', ...
    'LineWidth', 1, 'CapSize', 0)
hold on
errorbar(0.5, iti_Mm, iti_Mm_sem, 'ko', 'MarkerFaceColor', 'k', ...
    'LineWidth', 1, 'CapSize', 0)
errorbar(0.75, iti_Mh, iti_Mh_sem, 'ro', 'MarkerFaceColor', 'r', ...
    'LineWidth', 1, 'CapSize', 0)
xlim([0 1])
ylim([-0.2, 0.2])
xticks(0.25:0.25:0.75)
xticklabels({'Low', 'Mix', 'High'})
xlabel('Reward block')
ylabel('Mean trial initiation time (s)')
set(gca, 'TickDir', 'out'); box off;
% axis square
ax10 = gca;
ax10.YRuler.TickLabelGapOffset = 1;
title('\rm Muscimol')

%delta iti histogram
nexttile
h = histogram(deltaiti, bins);
h.FaceColor = [0.65 0.65 0.65];
h.EdgeColor = [0.65 0.65 0.65];
% text(-0.4, 2, strcat('p = ', num2str(pratio_iti)))
xlim([-0.5 0.5])
ylim([0 3])
xline(0, 'k--', 'LineWidth', 1)
ylabel('N rats')    
xlabel({'\Delta initiation time ratio'; '(muscimol - control)'})
% axis square
set(gca, 'TickDir', 'out'); box off;
ax11 = gca;
ax11.YRuler.TickLabelGapOffset = 1;


set(gcf,'renderer','painter')
set(gcf, 'units', 'centimeters', 'position', [10 10 16 4])
