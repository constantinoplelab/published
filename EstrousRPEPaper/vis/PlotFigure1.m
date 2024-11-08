function PlotFigure1(datadir, codedir)
%PlotFigure1 - Plots Figure 1. Must be run after ProcessData_Figure1.
% INPUTS:
%   datadir - Local directory where ProcessData_Figure1.mat was saved after running ProcessData_Figure1
%   codedir - Local directory where code (e.g., published/EstrousRPE/) was saved

%% Set up paths
s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codedir, 'IgnoreCase', ispc);
if ~onPath % only add code dir to path if it already isn't
    addpath(genpath(codedir))
end

%% Load data
load([datadir, 'ProcessData_Figure1'],...
    'f_ratlist', 'm_ratlist', 'ITIbyBlock', 'beh_sens_staged',...
    'newfrats', 'beh_sens_males', 'newmrats', 'Latencydynamics',...
    'delta_rats', 'E2_rats');

%% Create variables
cycle = {'Proestrus', 'Diestrus'};
group_colors = {'#E87003'; '#7358C6'; '#9E9D9D'}; %dark orange, purple, black
nback = 7;
frats = length(f_ratlist);
mrats = length(m_ratlist);

%% PLOT %%
figure;
set(gcf, 'Color', [1 1 1], 'Units', 'Inches', 'Position', [0, 0, 17, 9],...
    'PaperUnits', 'Inches', 'PaperSize', [9, 5], 'renderer','painters')   

% Means are used to summarize individual rat data
% Medians are used to compare across groups (unless noted otherwise)
% Statistical tests match means/medians (parametric/non-parametric)
% Detrended initiation times are compared across blocks and are regressed 
% against previous trial reward volumes
%--------------------------------------------------------------------------
%1c. Mean trial initiation times for low and high blocks for example rat
%--------------------------------------------------------------------------
nexttile
ratnum = find(strcmp('G027', f_ratlist));
pvals = [ITIbyBlock.Female.ranksum_pval_det ITIbyBlock.Male.ranksum_pval_det];
lo = ITIbyBlock.Female.low_ITI_det{ratnum};
mix = ITIbyBlock.Female.mixed_ITI_det{ratnum};
hi = ITIbyBlock.Female.high_ITI_det{ratnum};
lo_err = ITIbyBlock.Female.low_ITI_det_err{ratnum};
mix_err = ITIbyBlock.Female.mixed_ITI_det_err{ratnum};
hi_err = ITIbyBlock.Female.high_ITI_det_err{ratnum};
plot(1, lo, '.k', 'MarkerSize', 30); hold on
plot(2, mix, '.k', 'MarkerSize', 30); hold on
plot(3, hi, '.k', 'MarkerSize', 30); hold on
errorbar(1, lo, lo_err*1.96, color='k', LineWidth=0.5, capsize=10); hold on    
errorbar(2, mix, mix_err*1.96, color='k', LineWidth=0.5, capsize=10); hold on    
errorbar(3, hi, hi_err*1.96, color='k', LineWidth=0.5, capsize=10); hold on    
ylim([-1.5 1.5])
yticks(-2:1:2)
xlim([0.5 3.5]);
xticks(1:3)
xticklabels({'Low', 'Mixed', 'High'})    
ylabel('Initiation (init.) time (s, detrended)')
grid off
set(gca, 'TickDir', 'out'); box off
axis square 
subtitle(['Rat G027 p=' num2str(pvals(ratnum))])
title('c')

%--------------------------------------------------------------------------        
%1d. Distribution of differences in trial initiation times for low and high blocks for all rats
%--------------------------------------------------------------------------
nexttile
deltas = [ITIbyBlock.Female.deltas_det ITIbyBlock.Male.deltas_det];
histogram(deltas,'BinWidth', 0.25,...
    'FaceColor', 'k', 'EdgeColor', 'none'); hold on
alpha(0.2)
grid off
set(gca, 'TickDir', 'out'); box off
axis square 
ylabel('# rats')
xlabel('Low - high trial initiation')
xlim([-4.7 4.7])
ylim([0 60])
yticks(0:20:60)
ylimits = ylim;
plot([0 0],[0 ylimits(2)], '--', 'LineWidth', 0.75, 'Color', [0 0 0])
xline(median(deltas, 'omitnan'), '--r', 'LineWidth', 0.75) 
pval = signrank(deltas);
subtitle(['N=' num2str(length(deltas))...
    ', sign rank diff from zero, p=' num2str(pval)])
sigrats = pvals(pvals<0.05 & deltas>0);
disp(['d ' num2str(length(sigrats))...
    ' rats out of ' num2str(mrats + frats)...
    ' had significantly higher initiation times in low blocks than high blocks'])
title('d')

%--------------------------------------------------------------------------
%1f. Regression of previous trial reward volumes on detrended trial initiation times for all rats
%--------------------------------------------------------------------------
nexttile
clear l
betas = [ITIbyBlock.Female.betas_detrended; ITIbyBlock.Male.betas_detrended];
l(1) = shadedErrorBar(1:nback, median(betas(:, 1:nback), 'omitnan'),...
    sem(betas(:, 1:nback))*1.96, 'lineProps', {'-', 'Color', 'k',...
    'LineWidth', 0.5}); hold on 
xlim([0.5 nback-1+0.5]) 
xticks(1:nback)
xticklabels(1:nback)
xlabel('Trials back')
ylabel('Regression coefficient')
yticks(-0.3:0.1:0)
ylim([-0.3 0.05])
yline(0, 'k--', 'linewidth', 0.75)
set(gca, 'TickDir', 'out'); box off;
axis square
arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), l)
title('f')

%--------------------------------------------------------------------------
%1g. Sensitivity of detrended trial initiation times to blocks for all female rats in 
% proestrus and diestrus stages and male rats
%--------------------------------------------------------------------------
nexttile
clear l
pro_low = cell2mat(ITIbyBlock.(cycle{1}).low_ITI_det)';
pro_high = cell2mat(ITIbyBlock.(cycle{1}).high_ITI_det)';
di_low = cell2mat(ITIbyBlock.(cycle{2}).low_ITI_det)';
di_high = cell2mat(ITIbyBlock.(cycle{2}).high_ITI_det)';
male_low = cell2mat(ITIbyBlock.('Male').low_ITI_det)';
male_high = cell2mat(ITIbyBlock.('Male').high_ITI_det)';
l(1) = shadedErrorBar([1 2],...
    median([pro_low pro_high], 'omitnan'),...
    sem([pro_low, pro_high]),...
    'lineProps', {'-', 'Color', group_colors{1}, 'LineWidth', 0.5}); hold on
l(1).mainLine.DisplayName = cycle{1};
l(2) = shadedErrorBar([1 2],...
    median([di_low, di_high], 'omitnan'),...
    sem([di_low, di_high]),...
    'lineProps', {'-', 'Color', group_colors{2}, 'LineWidth', 0.5}); hold on
l(2).mainLine.DisplayName = cycle{2};
l(3) = shadedErrorBar([1 2],...
    median([male_low, male_high], 'omitnan'),...
    sem([male_low, male_high]),...
    'lineProps', {'-', 'Color', 'k', 'LineWidth', 0.5}); hold on
l(3).mainLine.DisplayName = 'Males';
legend('Location','best')
arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), l)
ylim([-0.9 0.8])
yticks(-1:0.5:1)
xlim([0.5 2.5])
xticks(1:2)
xticklabels({'Low', 'High'})
ylabel('Init. time (s, detrended)')
grid off; set(gca, 'TickDir', 'out'); box off; axis square
title('g')

%stats
low_x = [pro_low; di_low; male_low];
high_x = [pro_high; di_high; male_high];
group = char(repmat('Proestrus', length(pro_low), 1),...
    repmat('Diestrus', length(di_low), 1),...
    repmat('Male', length(male_low), 1));
pval_kruskal_low = kruskalwallis(low_x, group, 'off');
pval_kruskal_high = kruskalwallis(high_x, group, 'off');
pval_low_dipro = signrank(pro_low, di_low);
effsize_low_dipro = effsize(pro_low, di_low);
pval_high_dipro = signrank(pro_high, di_high);
effsize_high_dipro = effsize(pro_high, di_high);
pval_low_dimale = ranksum(di_low, male_low);
pval_high_dimale = ranksum(di_high, male_high);
pval_low_promale = ranksum(pro_low, male_low);
pval_high_promale = ranksum(pro_high, male_high);
subtitle([{'Kruskal wallis:'}; {['low p=' num2str(pval_kruskal_low)]}; {['high p=' num2str(pval_kruskal_high)]};...
    {'di vs pro:'}; {['sign rank low p=' num2str(pval_low_dipro) ', eff: ' num2str(effsize_low_dipro)]};...
        {['sign rank high p=' num2str(pval_high_dipro) ', eff: ' num2str(effsize_high_dipro)]};...
    {'di vs male:'}; {['rank sum low p=' num2str(pval_low_dimale)]};...
        {['rank sum high p=' num2str(pval_high_dimale)]};...
    {'pro vs male:'}; {['rank sum low p=' num2str(pval_low_promale)]};...
        {['rank sum high p=' num2str(pval_high_promale)]}])

%--------------------------------------------------------------------------
%1h. Plot as low - high, individual points and means
%--------------------------------------------------------------------------
nexttile
pro_deltas = ITIbyBlock.(cycle{1}).deltas_det;
di_deltas = ITIbyBlock.(cycle{2}).deltas_det;
male_deltas = ITIbyBlock.Male.deltas_det;
for rat = 1:length(f_ratlist)
    plot([1 2], [pro_deltas(rat) di_deltas(rat)],...
        '-', color = [0 0 0 0.1], linewidth = 0.5); hold on
end
for rat = 1:length(m_ratlist)
    plot(3+randn(1)*0.05, male_deltas(rat), '.',...
        color = [0.7 0.7 0.7], markersize = 10);
end
plot(1:3, [median(pro_deltas, 'omitnan')...
    median(di_deltas, 'omitnan')...
    median(male_deltas, 'omitnan')], '.k', markersize = 30); hold on
errorbar(1, median(pro_deltas, 'omitnan'), sem(pro_deltas),...
    color='k', LineWidth=0.5, capsize=10); hold on    
errorbar(2, median(di_deltas, 'omitnan'), sem(di_deltas),...
    color='k', LineWidth=0.5, capsize=10); hold on   
errorbar(3, median(male_deltas, 'omitnan'), sem(male_deltas),...
    color='k', LineWidth=0.5, capsize=10); hold on   
xlim([0.5 3.5]);
xticks(1:1:3)
xticklabels({'Proestrus', 'Diestrus', 'Males'})
ylabel('\Delta init. time')
% ylim([0 22])
yticks(-2:2:6)
grid off; set(gca, 'TickDir', 'out'); box off; axis square
pval_sr_ProDi_delta = signrank(pro_deltas, di_deltas);
pval_rs_ProMale_delta = ranksum(pro_deltas, male_deltas);
pval_rs_DiMale_delta = ranksum(di_deltas, male_deltas);
effsize_dipro = effsize(pro_deltas, di_deltas);
disp(['h Pro vs Di sign rank p=' num2str(pval_sr_ProDi_delta)])
disp(['h Pro vs Male rank sum p=' num2str(pval_rs_ProMale_delta)])
disp(['h Di vs Male rank sum p=' num2str(pval_rs_DiMale_delta)])
disp(['h Pro vs Di effect size=' num2str(effsize_dipro)])
title('h')

%--------------------------------------------------------------------------
%1i. Detrended initiation time dynamics across block transitions
%--------------------------------------------------------------------------
nexttile
twin = 30;
xvec = (-twin:1:twin);
set(gcf,'color','w','renderer','painters');
lowcolors = {'-b', '-k'};
highcolors = {'-r', '-k'};
hightolow_stages = cell(length(cycle), 1);
lowtohigh_stages = cell(length(cycle), 1);
for e = 1:length(cycle)

    ltom = Latencydynamics.(cycle{e}).ltom;
    htom = Latencydynamics.(cycle{e}).htom;
    mtol = Latencydynamics.(cycle{e}).mtol;
    mtoh = Latencydynamics.(cycle{e}).mtoh;

    % Plot as value transitions
    %group high to low value transitions
    hightolow = [htom; mtol]; %high to mix, mix to low
    hightolow_median = median(hightolow, 'omitnan');
    hightolow_stages{e} = hightolow;

    %group low to high value transitions
    %low to mix, mix to high
    lowtohigh = [ltom; mtoh];
    lowtohigh_median = median(lowtohigh, 'omitnan');
    lowtohigh_stages{e} = lowtohigh;
    clear l
    l(1) = shadedErrorBar(xvec, hightolow_median, ...
        std(hightolow,...
        'omitnan')./sqrt(length(hightolow(:,1))/2),...
        'lineprops', lowcolors{e});
    l(2) = shadedErrorBar(xvec, lowtohigh_median, ...
        std(lowtohigh,...
        'omitnan')./sqrt(length(lowtohigh(:,1))/2),...
        'lineprops', highcolors{e});
    arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), l)
    hold on
    set(gca, 'TickDir', 'out'); box off
    xlabel('Trials from block switch');
    ylabel('Init. time (s, detrended)'); %removed z-score
    
    xlim([-twin/3 twin])
    yticks(-2:1:2)
    axis square
    set(gcf, 'Color', [1 1 1]);

end
thisrange = ylim;
line([0 0], [thisrange(1) thisrange(2)], 'Color', [0 0 0], 'LineStyle', '--');
title('i High to low value/Low to high value');

%Statistics
%High value to low value
range = ylim;
hightolow_pro = hightolow_stages{1};
hightolow_di = hightolow_stages{2};
pvals_hightolow = NaN(size(hightolow_pro, 2), 1);
for t = 1:size(hightolow_pro, 2) %loop over trials from transition
    pvals_hightolow(t) = signrank(hightolow_pro(:, t), hightolow_di(:, t));
    if pvals_hightolow(t) < 0.05 %add significance stars
        text(xvec(t), range(2)+0.2, '*', color = 'b')
    end
end

%Low value to high value
lowtohigh_pro = lowtohigh_stages{1};
lowtohigh_di = lowtohigh_stages{2};
pvals_lowtohigh = NaN(size(lowtohigh_pro, 2), 1);
for t = 1:size(lowtohigh_pro, 2) %loop over trials from transition
    pvals_lowtohigh(t) = signrank(lowtohigh_pro(:, t), lowtohigh_di(:, t));
    if pvals_lowtohigh(t) < 0.05 %add significance stars
        text(xvec(t), range(1)+0.2, '*', color = 'r')
    end
end

%--------------------------------------------------------------------------
%1j. Regression weights of trial initiation times against rewards on previous trials
% for female rats in proestrus and diestrus stages and males. Excludes rats
% that don't have greater initiation times in low compared to high blocks.
%--------------------------------------------------------------------------
nexttile
clear l
pro_betas = ITIbyBlock.(cycle{1}).betas_detrended(beh_sens_staged, :);
di_betas = ITIbyBlock.(cycle{2}).betas_detrended(beh_sens_staged, :);
male_betas = ITIbyBlock.('Male').betas_detrended(beh_sens_males, :);

l(1)=shadedErrorBar(1:nback, mean(pro_betas(:, 1:nback), 'omitnan'),...
    sem(pro_betas(:, 1:nback))*1.96, 'lineProps', {'-', 'Color',...
    group_colors{1}, 'LineWidth', 0.5}); hold on
l(1).mainLine.DisplayName = cycle{1};
l(2)=shadedErrorBar(1:nback, mean(di_betas(:, 1:nback), 'omitnan'),...
    sem(di_betas(:, 1:nback))*1.96, 'lineProps', {'-', 'Color',...
    group_colors{2}, 'LineWidth', 0.5}); hold on
l(2).mainLine.DisplayName = cycle{2};
l(3)=shadedErrorBar(1:nback, mean(male_betas(:, 1:nback), 'omitnan'),...
    sem(male_betas(:, 1:nback))*1.96, 'lineProps', {'-', 'Color',...
    'k', 'LineWidth', 0.5}); hold on
l(3).mainLine.DisplayName = 'Males';
legend('Location','best')
xlim([0.5 6.5])
ylim([-0.65 0.05])
yticks(-0.8:0.2:0.8)
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
     % regresspvals(prevrew, 1) = kruskalwallis([pro_betas(:, prevrew);...
     %    di_betas(:, prevrew); male_betas(:, prevrew)], [ones(newfrats, 1);...
     %    ones(newfrats, 1)*2; ones(newmrats, 1)*3], 'off');     
    regresspvals(prevrew, 1) = anova1([pro_betas(:, prevrew);...
        di_betas(:, prevrew); male_betas(:, prevrew)], [ones(newfrats, 1);...
        ones(newfrats, 1)*2; ones(newmrats, 1)*3], 'off');
    if regresspvals(prevrew, 1)<0.05
        text(prevrew, 0.02, '*', 'FontSize',14)
    end    
end
effsize_t_1 = effsize(pro_betas(:, 1), di_betas(:, 1));
disp('j anova p-values: ')
disp(regresspvals)
disp(['j effect size p vs d t-1: ' num2str(effsize_t_1)])
title('j')

%--------------------------------------------------------------------------
%1k.  Coefficients for one trial back for proestrus vs diestrus stages 
%--------------------------------------------------------------------------
nexttile
pro_betas = ITIbyBlock.(cycle{1}).betas_detrended(beh_sens_staged, :);
di_betas = ITIbyBlock.(cycle{2}).betas_detrended(beh_sens_staged, :);
stageeffect_beta_Rt_1 = di_betas(:, 1) - pro_betas(:, 1);
these_binedges = [-2.05:0.1:-0.05 0.05:0.1:7.05];
histogram(stageeffect_beta_Rt_1, 'FaceColor', 'k', binedges = these_binedges); hold on
xline(0, '--k', linewidth=0.75)
xline(median(stageeffect_beta_Rt_1, 'omitnan'), '--r', 'LineWidth', 0.75)
xlim([-1.5 1.5]);
xticks(-2:1:2)
xlabel({'Diestrus - proestrus beta Rt-1'})
ylabel('# rats')
ylim([0 26])
yticks(0:10:25)
set(gca, 'TickDir', 'out'); box off
pval_sr_ProDi_delta = signrank(di_betas(:, 1), pro_betas(:, 1));
title(['k Pro vs Di beta Rt-1 sign rank p=' num2str(pval_sr_ProDi_delta)])
grid off; set(gca, 'TickDir', 'out'); box off; axis square

%--------------------------------------------------------------------------
%1l.  Coefficients for one trial back for proestrus and diestrus stages as
%a function of hormonal modulation (proestrus - diestrus delta initiation times)
%--------------------------------------------------------------------------
nexttile
hormonal_modulation = pro_deltas(beh_sens_staged) - di_deltas(beh_sens_staged);
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

%--------------------------------------------------------------------------
%1m. Hormonal modulation as a function of changes in estradiol
%--------------------------------------------------------------------------
nexttile
x = E2_rats(:, 1) - E2_rats(:, 2); 
y = delta_rats(:, 1) - delta_rats(:, 2);
p = polyfit(x,y,1); 
yi = polyval(p,x) ; 
plot(x, y, '.', color='#9E9D9D', markersize = 15); hold on
plot(x,yi,'--k',linewidth=0.5); hold on
xlabel('Proestrus - diestrus estradiol (pg/ml)')
ylabel('Proestrus - diestrus \Delta init. time')
lsline
[R,P] = corrcoef(x, y); 
title(['m R=' num2str(R(2,1)) ', p=' num2str(P(2,1))])
ylim([-1.7 3.5])
yticks(-2:1:4)
xticks(15:20:55)
grid off; axis square; set(gca, 'TickDir', 'out'); box off

end