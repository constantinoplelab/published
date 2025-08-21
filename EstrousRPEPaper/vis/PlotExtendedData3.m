function PlotExtendedData3(datadir, codedir)
%PlotExtendedData1 - Plots Extended Data 3. 
% INPUTS:
%   datadir - Local directory where ProcessData_ExtendedData3.mat was saved after running ProcessData_ExtendedData3
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
load([datadir 'ProcessData_ExtendedData3'], 'ITIbyBlockPhotometry',...
    'ITIbyBlockSerum','NAcc_ratlist', 'Serum_ratlist');

%% Create variables
cycle = {'Proestrus', 'Diestrus'};
cyclecolors = {'#E87003'; '#7358C6'}; %dark orange, purple
nback = 7;

%% PLOT %%
figure;
set(gcf, 'Color', [1 1 1], 'Units', 'Inches', 'Position', [0, 0, 17, 9],...
    'PaperUnits', 'Inches', 'PaperSize', [9, 5], 'renderer','painters')

%--------------------------------------------------------------------------
%ED3a. Regression of previous trial reward volumes on trial initiation times
%for serum estradiol rats
%--------------------------------------------------------------------------
nexttile
clear l
betas = ITIbyBlockSerum.Female.betas_detrended;
l(1) = shadedErrorBar(1:nback, median(betas(:, 1:nback), 'omitnan'),...
    sem(betas(:, 1:nback))*1.96, 'lineProps', {'-', 'Color', 'k',...
    'LineWidth', 0.5}); hold on
xlim([0.5 nback-1+0.5]) 
xticks(1:nback)
xticklabels(1:nback)
xlabel('Trials back')
ylabel('Regression coefficient')
yticks(-1:0.1:0)
ylim([-0.4 0.05])
yline(0, 'k--', 'linewidth', 0.75)
set(gca, 'TickDir', 'out'); box off;
axis square
arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), l)
title(['a serum N=' num2str(length(Serum_ratlist))])

%--------------------------------------------------------------------------
% ED3b. Sensitivity of trial initiation times to blocks for all female rats in 
% proestrus and diestrus stages and male rats
%--------------------------------------------------------------------------
nexttile
clear l
pro_low = cell2mat(ITIbyBlockSerum.(cycle{1}).low_ITI_det)';
pro_high = cell2mat(ITIbyBlockSerum.(cycle{1}).high_ITI_det)';
di_low = cell2mat(ITIbyBlockSerum.(cycle{2}).low_ITI_det)';
di_high = cell2mat(ITIbyBlockSerum.(cycle{2}).high_ITI_det)';
l(1) = shadedErrorBar([1 2],...
    median([pro_low pro_high], 'omitnan'),...
    sem([pro_low, pro_high]),...
    'lineProps', {'-', 'Color', cyclecolors{1}, 'LineWidth', 0.5}); hold on
l(1).mainLine.DisplayName = cycle{1};
l(2) = shadedErrorBar([1 2],...
    median([di_low, di_high], 'omitnan'),...
    sem([di_low, di_high]),...
    'lineProps', {'-', 'Color', cyclecolors{2}, 'LineWidth', 0.5}); hold on
l(2).mainLine.DisplayName = cycle{2};
legend('Location','northeast')
arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), l)
xlim([0.5 2.5])
xticks(1:2)
xticklabels({'Low', 'High'})
ylabel('Initiation time (s)')
grid off; set(gca, 'TickDir', 'out'); box off; axis square
title('b serum')

%stats
pval_low_dipro = signrank(pro_low, di_low);
pval_high_dipro = signrank(pro_high, di_high);
subtitle(['Sign rank: low p=' num2str(pval_low_dipro) ', high p=' num2str(pval_high_dipro)])
effectsize_estrous_low = (mean(pro_low) - mean(di_low))./...
    std([pro_low di_low]); %for P vs D
disp(['P vs D low block effect size = ' num2str(effectsize_estrous_low)]) 
effectsize_estrous_high = (mean(pro_high) - mean(di_high))./...
    std([pro_high di_high]); %for P vs D
disp(['P vs D high block effect size = ' num2str(effectsize_estrous_high)]) 

%--------------------------------------------------------------------------
% Serum rats: Regression weights of trial initiation times against rewards on previous trials
% for female rats in proestrus and diestrus stages.
%--------------------------------------------------------------------------
nexttile
clear l
pro_betas = ITIbyBlockSerum.(cycle{1}).betas_detrended;
di_betas = ITIbyBlockSerum.(cycle{2}).betas_detrended;
l(1)=shadedErrorBar(1:nback, median(pro_betas(:, 1:nback), 'omitnan'),...
    sem(pro_betas(:, 1:nback))*1.96, 'lineProps', {'-', 'Color',...
    cyclecolors{1}, 'LineWidth', 0.5}); hold on
l(1).mainLine.DisplayName = cycle{1};
l(2)=shadedErrorBar(1:nback, median(di_betas(:, 1:nback), 'omitnan'),...
    sem(di_betas(:, 1:nback))*1.96, 'lineProps', {'-', 'Color',...
    cyclecolors{2}, 'LineWidth', 0.5}); hold on
l(2).mainLine.DisplayName = cycle{2};
legend('Location','best')
xlim([0.5 6.5])
ylim([-0.5 0.05])
yticks(-1:0.2:0.8)
xticks(1:7)
xticklabels(1:7)
yline(0, '--', color = 'k', linewidth=0.75, HandleVisibility='off')
xlabel('Trials Back')
ylabel('Regression coefficient')
grid off
set(gca, 'TickDir', 'out'); box off
axis square
grid off; set(gca, 'TickDir', 'out'); box off; axis square
arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), l)
regresspvals = NaN(nback, 1);
for prevrew = 1:nback  
     regresspvals(prevrew, 1) = signrank(pro_betas(:, prevrew),...
         di_betas(:, prevrew));     
    if regresspvals(prevrew, 1)<0.05
        text(prevrew, 0.02, '*', 'FontSize',14)
    end    
end
disp(regresspvals)
title('serum')

%--------------------------------------------------------------------------
% ED3c. Coefficients for one trial back for proestrus vs diestrus stages
%histogram
%--------------------------------------------------------------------------
nexttile
stageeffect_beta_Rt_1 = di_betas(:, 1) - pro_betas(:, 1);
these_binedges = [-0.56:0.08:-0.04 0.04:0.08:0.56];
histogram(stageeffect_beta_Rt_1, 'FaceColor', 'k', binedges = these_binedges); hold on
xline(0, '--k', linewidth=0.75)
xline(median(stageeffect_beta_Rt_1, 'omitnan'), '-k', linewidth=0.75)
xlim([-0.4 0.4]);
xticks(-0.4:0.2:0.4);
xlabel({'Diestrus - proestrus beta Rt-1'})
ylabel('# rats')
ylim([0 5.3])
yticks(1:5)
set(gca, 'TickDir', 'out'); box off
pval_sr_ProDi_delta = signrank(di_betas(:, 1), pro_betas(:, 1));
title('c serum')
subtitle(['Pro vs Di beta Rt-1 sign rank p=' num2str(pval_sr_ProDi_delta)])
grid off; set(gca, 'TickDir', 'out'); box off; axis square

%%%%%%%%%%%%%%%%%%%%%%%%%% PHOTOMETRY RATS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% ED3e. Regression of previous trial reward volumes on trial initiation times for photometry rats
%--------------------------------------------------------------------------
nexttile
clear l
betas = ITIbyBlockPhotometry.Female.betas_detrended;
l(1) = shadedErrorBar(1:nback, median(betas(:, 1:nback), 'omitnan'),...
    sem(betas(:, 1:nback))*1.96, 'lineProps', {'-', 'Color', 'k',...
    'LineWidth', 0.5}); hold on
xlim([0.5 nback-1+0.5]) 
xticks(1:nback)
xticklabels(1:nback)
xlabel('Trials back')
ylabel('Regression coefficient')
yticks(-1:0.2:0)
ylim([-1 0.05])
yline(0, 'k--', 'linewidth', 0.75)
set(gca, 'TickDir', 'out'); box off;
axis square
arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), l)
title(['e photometry N=' num2str(length(NAcc_ratlist))])

%--------------------------------------------------------------------------
% ED3f. Sensitivity of trial initiation times to blocks for all female rats in 
% proestrus and diestrus stages and male rats
%--------------------------------------------------------------------------
nexttile
clear l
pro_low = cell2mat(ITIbyBlockPhotometry.(cycle{1}).low_ITI_det)';
pro_high = cell2mat(ITIbyBlockPhotometry.(cycle{1}).high_ITI_det)';
di_low = cell2mat(ITIbyBlockPhotometry.(cycle{2}).low_ITI_det)';
di_high = cell2mat(ITIbyBlockPhotometry.(cycle{2}).high_ITI_det)';
l(1) = shadedErrorBar([1 2],...
    median([pro_low pro_high], 'omitnan'),...
    sem([pro_low, pro_high]),...
    'lineProps', {'-', 'Color', cyclecolors{1}, 'LineWidth', 0.5}); hold on
l(1).mainLine.DisplayName = cycle{1};
l(2) = shadedErrorBar([1 2],...
    median([di_low, di_high], 'omitnan'),...
    sem([di_low, di_high]),...
    'lineProps', {'-', 'Color', cyclecolors{2}, 'LineWidth', 0.5}); hold on
l(2).mainLine.DisplayName = cycle{2};
legend('Location','northeast')
arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), l)
ylim([-1.5 1.6])
yticks(-2:1:2)
xlim([0.5 2.5])
xticks(1:2)
xticklabels({'Low', 'High'})
ylabel('Init. time (s, detrended)')
grid off; set(gca, 'TickDir', 'out'); box off; axis square
title('f photometry')

%stats
pval_low_dipro = signrank(pro_low, di_low);
pval_high_dipro = signrank(pro_high, di_high);
subtitle(['Sign rank: low p=' num2str(pval_low_dipro) ', high p=' num2str(pval_high_dipro)])
effectsize_estrous_low = (mean(pro_low) - mean(di_low))./...
    std([pro_low di_low]); %for P vs D
disp(['P vs D low block effect size = ' num2str(effectsize_estrous_low)]) 
effectsize_estrous_high = (mean(pro_high) - mean(di_high))./...
    std([pro_high di_high]); %for P vs D
disp(['P vs D high block effect size = ' num2str(effectsize_estrous_high)]) 

%--------------------------------------------------------------------------
% Regression weights of trial initiation times against rewards on previous trials
% for female rats in proestrus and diestrus stages. 
%--------------------------------------------------------------------------
nexttile
clear l
pro_betas = ITIbyBlockPhotometry.(cycle{1}).betas_detrended;
di_betas = ITIbyBlockPhotometry.(cycle{2}).betas_detrended;
l(1)=shadedErrorBar(1:nback, median(pro_betas(:, 1:nback), 'omitnan'),...
    sem(pro_betas(:, 1:nback))*1.96, 'lineProps', {'-', 'Color',...
    cyclecolors{1}, 'LineWidth', 0.5}); hold on
l(1).mainLine.DisplayName = cycle{1};
l(2)=shadedErrorBar(1:nback, median(di_betas(:, 1:nback), 'omitnan'),...
    sem(di_betas(:, 1:nback))*1.96, 'lineProps', {'-', 'Color',...
    cyclecolors{2}, 'LineWidth', 0.5}); hold on
l(2).mainLine.DisplayName = cycle{2};
legend('Location','best')
xlim([0.5 6.5])
ylim([-1 0.05])
yticks(-1:0.2:0.8)
xticks(1:7)
xticklabels(1:7)
yline(0, '--', color = 'k', linewidth=0.75, HandleVisibility='off')
xlabel('Trials Back')
ylabel('Regression coefficient')
grid off
set(gca, 'TickDir', 'out'); box off
axis square
grid off; set(gca, 'TickDir', 'out'); box off; axis square
arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), l)
regresspvals = NaN(nback, 1);
for prevrew = 1:nback  
     regresspvals(prevrew, 1) = signrank(pro_betas(:, prevrew),...
         di_betas(:, prevrew));     
    if regresspvals(prevrew, 1)<0.05
        text(prevrew, 0.02, '*', 'FontSize',14)
    end    
end
disp(regresspvals)
title('photometry')

%--------------------------------------------------------------------------
%ED3f. Coefficients for one trial back for proestrus vs diestrus stages
%histogram
%--------------------------------------------------------------------------
nexttile
stageeffect_beta_Rt_1 = di_betas(:, 1) - pro_betas(:, 1);
these_binedges = [-0.56:0.08:-0.04 0.04:0.08:0.56];
histogram(stageeffect_beta_Rt_1, 'FaceColor', 'k', binedges = these_binedges); hold on
xline(0, '--k', linewidth=0.75)
xline(median(stageeffect_beta_Rt_1, 'omitnan'), '-k', linewidth=0.75)
xlim([-0.5 0.5]);
xticks(-0.5:0.5:0.5);
xlabel({'Diestrus - proestrus beta Rt-1'})
ylabel('# rats')
ylim([0 3.2])
yticks(1:4)
set(gca, 'TickDir', 'out'); box off
pval_sr_ProDi_delta = signrank(di_betas(:, 1), pro_betas(:, 1));
title('f photometry')
subtitle(['Pro vs Di beta Rt-1 sign rank p=' num2str(pval_sr_ProDi_delta)])
grid off; set(gca, 'TickDir', 'out'); box off; axis square


end