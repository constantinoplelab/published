function PlotExtendedData1(datadir, codedir) 
%PlotExtendedData1 - Plots Extended Data 1. 
% INPUTS:
%   datadir - Local directory where ProcessData_ExtendedData1.mat was saved after running ProcessData_ExtendedData1
%   codedir - Local directory where code (e.g., published/EstrousRPE/) was saved

%% Set up paths
s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codedir, 'IgnoreCase', ispc);
if ~onPath % only add code dir to path if it isn't already
    addpath(genpath(codedir))
end

%% Load data
load([datadir 'ProcessData_ExtendedData1'], 'f_ratlist',...
    'Stages', 'Concentration', 'SerumTable','ITIbyBlock',...
    'beh_sens_staged');

%% PLOT %%
figure;
set(gcf, 'Color', [1 1 1], 'Units', 'Inches', 'Position', [0, 0, 17, 9],...
    'PaperUnits', 'Inches', 'PaperSize', [9, 5], 'renderer','painters')

%--------------------------------------------------------------------------
%ED1c. Estradiol serum levels, with rats as individual points
%--------------------------------------------------------------------------
nexttile
errorbar(1:4, [mean(Concentration.D, 'omitnan') median(Concentration.PE, 'omitnan')...
    median(Concentration.PL, 'omitnan') median(Concentration.E, 'omitnan')],...
    [std(Concentration.D, 'omitnan')/sqrt(sum(~isnan(Concentration.D)))...
    std(Concentration.PE, 'omitnan')/sqrt(sum(~isnan(Concentration.PE)))...
    std(Concentration.PL, 'omitnan')/sqrt(sum(~isnan(Concentration.PL)))...
    std(Concentration.E, 'omitnan')/sqrt(sum(~isnan(Concentration.E)))], '.k', markersize=30,...
    linewidth=0.5, capsize=10); hold on
for d=1:length(Concentration.D)
    plot(0.8+randn(1)*0.03, Concentration.D(d), '.', color='#9E9D9D', markersize=15); hold on
end
for pe=1:length(Concentration.PE)
    plot(1.8+randn(1)*0.03, Concentration.PE(pe), '.', color='#9E9D9D', markersize=15); hold on
end
for pl=1:length(Concentration.PL)
    plot(2.8+randn(1)*0.03, Concentration.PL(pl), '.', color='#9E9D9D', markersize=15); hold on    
end
for e=1:length(Concentration.E)
    plot(3.8+randn(1)*0.03, Concentration.E(e), '.', color='#9E9D9D', markersize=15); hold on    
end
xticks(1:4)
xlim([0.5 4.5])
ylim([0 12.5])
yticks(0:4:12)
xticklabels({'Diestrus', 'Proestrus (morning)',...
    'Proestrus (evening)', 'Estrus'})
ylabel('Estradiol (pg/ml)')
%statistical tests comparing expression over cycle
pval_kw = kruskalwallis([Concentration.D, [Concentration.PE; NaN],...
    [Concentration.PL; NaN], Concentration.E],[],'off');
pval_PLvD = ranksum(Concentration.D, Concentration.PL);
pval_PEvD = ranksum(Concentration.D, Concentration.PE);
pval_EvD = ranksum(Concentration.D, Concentration.E);
pval_PLvE = ranksum(Concentration.PL, Concentration.E);
pval_PEvE = ranksum(Concentration.PE, Concentration.E);
pval_PEvPL = ranksum(Concentration.PE, Concentration.PL);
effectsize_PvD = effsize(Concentration.D, Concentration.PL);
title(['c kruskal wallis p=' num2str(pval_kw)])
subtitle({['early proestrus vs late proestrus rank sum p=' num2str(pval_PEvPL)],...
    ['estrus vs early proestrus rank sum p=' num2str(pval_PEvE)],...
    ['estrus vs late proestrus rank sum p=' num2str(pval_PLvE)],...
    ['diestrus vs estrus rank sum p=' num2str(pval_EvD)],...
    ['diestrus vs early proestrus rank sum p=' num2str(pval_PEvD)],...
    ['diestrus vs late proestrus rank sum p=' num2str(pval_PLvD)],...
    ['diestrus vs late proestrus effect size=' num2str(effectsize_PvD)]})
set(gca, 'TickDir', 'out'); box off; grid off

%--------------------------------------------------------------------------
%ED1d. Median serum estradiol by stage of training rats that
% got blood draws
%--------------------------------------------------------------------------
%calculate avg per stage per rat
ratlist = unique(SerumTable.Rat); %training rats that got blood draws
estradiol_by_rat = NaN(length(ratlist), 2);
for rat = 1:length(ratlist)
    estradiol_by_rat(rat, 1) = median(SerumTable.EstradiolConc(strcmp(SerumTable.Stage, 'Proestrus')...
        & strcmp(SerumTable.Rat, ratlist{rat}) & SerumTable.Volume==50), 'omitnan');
    estradiol_by_rat(rat, 2) = median(SerumTable.EstradiolConc(strcmp(SerumTable.Stage, 'Diestrus')...
        & strcmp(SerumTable.Rat, ratlist{rat}) & SerumTable.Volume==50), 'omitnan');
end
nexttile
for rat = 1:length(ratlist)
    plot([1 2], estradiol_by_rat(rat, :), '-', color = [0 0 0 0.25],...
        linewidth=1); hold on
end
plot(1, median(estradiol_by_rat(:, 1)), '.k', markersize=30); hold on   
plot(2, median(estradiol_by_rat(:, 2)), '.k', markersize=30); hold on   
errorbar(1, median(estradiol_by_rat(:, 1)), sem(estradiol_by_rat(:, 1)),...
    'k', capsize=10, linewidth=0.5); hold on   
errorbar(2, median(estradiol_by_rat(:, 2)), sem(estradiol_by_rat(:, 2)),...
    'k', capsize=10, linewidth=0.5); hold on   
xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Proestrus', 'Diestrus'})
xlabel('Training session cycle stage')
ylim([20 220])
ylabel('Estradiol (pg/ml)')
pval = signrank(estradiol_by_rat(:, 1), estradiol_by_rat(:, 2));
effectsize = effsize(estradiol_by_rat(:, 1), estradiol_by_rat(:, 2));
title(['d sign rank p=' num2str(pval) ', effect size=' num2str(effectsize)])
grid off; axis square; set(gca, 'TickDir', 'out'); box off

%--------------------------------------------------------------------------
%ED1e. Sensitivity of detrended trial initiation times to blocks for
% female rats in each stage of the estrous cycle
%--------------------------------------------------------------------------
nexttile
frats = length(f_ratlist);
deltas_over_rats = NaN(frats, length(Stages));
for e=1:length(Stages)
    deltas = ITIbyBlock.(Stages{e}).deltas_det';
    errorbar(e, median(deltas, 'omitnan'), sem(deltas), 'k', capsize=10); hold on
    plot(e, median(deltas, 'omitnan'), '.k', markersize=30); hold on
    deltas_over_rats(:,e) = deltas;
end
xlim([0.5 4.5])
xticks(1:4)
xticklabels(Stages)
yticks(0.6:0.2:2)
ylabel('\Delta init. time (detrended)')
grid off; set(gca, 'TickDir', 'out'); box off; axis square
p = kruskalwallis(deltas_over_rats(:),...
    [ones(frats,1); 2*ones(frats,1);...
    3*ones(frats,1); 4*ones(frats,1)],...
    'off');
title(['e kruskal wallis p = ' num2str(p)])
%posthoc tests
MvsDpval = signrank(deltas_over_rats(:,1), deltas_over_rats(:,2));
PvsDpval = signrank(deltas_over_rats(:,2), deltas_over_rats(:,3));
PvsEpval = signrank(deltas_over_rats(:,3), deltas_over_rats(:,4));
PvsMpval = signrank(deltas_over_rats(:,1), deltas_over_rats(:,3));
EvsMpval = signrank(deltas_over_rats(:,1), deltas_over_rats(:,4));
EvsDpval = signrank(deltas_over_rats(:,2), deltas_over_rats(:,4));
disp(['d posthoc: sign rank P vs. E p=' num2str(PvsEpval),...
    ', sign rank E vs. M p=' num2str(EvsMpval),...
    ', sign rank M vs. D p=' num2str(MvsDpval),...
    ', sign rank P vs. M p=' num2str(PvsMpval),...
    ', sign rank P vs. D p=' num2str(PvsDpval),...
    ', sign rank E vs. D p=' num2str(EvsDpval)])

%--------------------------------------------------------------------------
%ED1f.  Coefficients for one trial back for proestrus and diestrus stages as
%a function of hormonal modulation (proestrus - diestrus delta initiation times)
%--------------------------------------------------------------------------
nexttile
cycle = {'Proestrus', 'Diestrus'};
pro_deltas = ITIbyBlock.(cycle{1}).deltas_det;
di_deltas = ITIbyBlock.(cycle{2}).deltas_det;
hormonal_modulation = pro_deltas(beh_sens_staged) - di_deltas(beh_sens_staged);
pro_betas = ITIbyBlock.(cycle{1}).betas_detrended(beh_sens_staged, :);
di_betas = ITIbyBlock.(cycle{2}).betas_detrended(beh_sens_staged, :);
stageeffect_beta_Rt_1 = di_betas(:, 1) - pro_betas(:, 1);
nan_idx = isnan(hormonal_modulation') | isnan(stageeffect_beta_Rt_1);
x = hormonal_modulation(~nan_idx);
y = stageeffect_beta_Rt_1(~nan_idx);
p = polyfit(x,y,1); 
yi = polyval(p,x); 
plot(x, y, '.', color='#9E9D9D', markersize = 15); hold on
plot(x,yi,'--k',linewidth=0.5)
xlabel('/Delta Initiation time (proestrus - diestrus)'); 
ylabel('Rt-1 (diestrus - proestrus)')
ylim([-1.2 1.5])
yticks(-1:1:2)
xlim([-2.5 4.5])
xticks(-4:2:6)
[R,P] = corrcoef(x', y);
title(['l R=' num2str(R(2,1)) ', p=' num2str(P(2,1))])
grid off; set(gca, 'TickDir', 'out'); box off; axis square


end