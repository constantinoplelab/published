function PlotFigure6(datadir, codedir)
%PlotFigure6 - Plots Figure 6. 
% INPUTS:
%   datadir - Local directory where ProcessData_Figure6.mat was saved after running ProcessData_Figure6
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
load([datadir, 'ProcessData_Figure6'], 'ratlist_shEsr1',...
    'TABLE_tet', 'TABLE_lenti', 'ITIbyBlock_shEsr1_screen',...
    'ITIbyBlock_shEsr1_noscreen', 'ratTrial', 'Ls', 'estradiol',...
    'RatConc');

%% PLOT %%
figure;
set(gcf, 'Color', [1 1 1], 'Units', 'Inches', 'Position', [0, 0, 17, 9],...
    'PaperUnits', 'Inches', 'PaperSize', [9, 5], 'renderer','painters')

%set generally used variables%% 
doxstatecolors = {'k', '#47C90A'};

%--------------------------------------------------------------------------
% 6c. Quantification of RNAscope for tet-shEsr1
%--------------------------------------------------------------------------
nexttile
%Controls
con_idx = contains(TABLE_tet.Group,'con') & contains(TABLE_tet.Probe,'Esr1');
con_perimage = TABLE_tet.Number_of_spots(con_idx)./...
    TABLE_tet.Number_of_images(con_idx);
con_no = con_perimage(con_perimage<(mean(con_perimage, 'omitnan')+2*std(con_perimage, 'omitnan')));
con_norm = con_no./(mean(con_no, 'omitnan')) * 100;
for pt = 1:length(con_norm)
    plot(1+randn(1)*0.05, con_norm(pt), '.', color = '#9E9D9D',...
        markersize=15); hold on
end
errorbar(1, median(con_norm, 'omitnan'), sem(con_norm),...
    'k', CapSize=10, linewidth=0.5)
plot(1, median(con_norm, 'omitnan'), '.k', markersize=30)
%shEsr1
exp_idx = contains(TABLE_tet.Group,'exp') & contains(TABLE_tet.Probe,'Esr1');
exp_perimage = TABLE_tet.Number_of_spots(exp_idx)./...
    TABLE_tet.Number_of_images(exp_idx);
exp_no = exp_perimage(exp_perimage<(mean(exp_perimage, 'omitnan')+2*std(exp_perimage, 'omitnan')));
exp_norm = exp_no./(mean(con_no, 'omitnan'))* 100;
for pt = 1:length(exp_norm)
    plot(2+randn(1)*0.05, exp_norm(pt), '.', color = '#9E9D9D',...
        markersize=15); hold on
end
errorbar(2, median(exp_norm, 'omitnan'), sem(exp_norm),...
    'k', CapSize=10, linewidth=0.5)
plot(2, median(exp_norm, 'omitnan'), '.k', markersize=30)
xlim([0.5 2.5])
xticks(1:2)
xticklabels({'Control', 'shEsr1'})
ylabel('% change in Esr1+ spots/image')
yline(100, '--k')
yticks(0:50:200)
ylim([0 200])
axis square; set(gca, 'TickDir', 'out'); box off; grid off
pval = ranksum(con_norm, exp_norm);
effectsize = (mean(con_norm, 'omitnan') - mean(exp_norm, 'omitnan'))...
    ./std([con_norm; exp_norm]);
title(['c rank sum p=' num2str(pval) ', effect size =' num2str(effectsize)])

%--------------------------------------------------------------------------
% 6d. Quantification of RNAscope for lenti-shEsr1
%--------------------------------------------------------------------------
nexttile
%Controls
con_idx = contains(TABLE_lenti.Group,'con') & contains(TABLE_lenti.Probe,'Esr1');
con_perimage = TABLE_lenti.Number_of_spots(con_idx)./...
    TABLE_lenti.Number_of_images(con_idx);
con_no = con_perimage(con_perimage<(mean(con_perimage, 'omitnan')+2*std(con_perimage, 'omitnan')));
con_norm = con_no./(mean(con_no, 'omitnan')) * 100;
for pt = 1:length(con_norm)
    plot(1+randn(1)*0.05, con_norm(pt), '.', color = '#9E9D9D',...
        markersize=15); hold on
end
errorbar(1, median(con_norm, 'omitnan'), sem(con_norm),...
    'k', CapSize=10, linewidth=0.5)
plot(1, median(con_norm, 'omitnan'), '.k', markersize=30)
%shEsr1
exp_idx = contains(TABLE_lenti.Group,'exp') & contains(TABLE_lenti.Probe,'Esr1');
exp_perimage = TABLE_lenti.Number_of_spots(exp_idx)./...
    TABLE_lenti.Number_of_images(exp_idx);
exp_no = exp_perimage(exp_perimage<(mean(exp_perimage, 'omitnan')+2*std(exp_perimage, 'omitnan')));
exp_norm = exp_no./(mean(con_no, 'omitnan')) * 100;
for pt = 1:length(exp_norm)
    plot(2+randn(1)*0.05, exp_norm(pt), '.', color = '#9E9D9D',...
        markersize=15); hold on
end
errorbar(2, median(exp_norm, 'omitnan'), sem(exp_norm),...
    'k', CapSize=10, linewidth=0.5)
plot(2, median(exp_norm, 'omitnan'), '.k', markersize=30)
xlim([0.5 2.5])
xticks(1:2)
xticklabels({'Control', 'shEsr1'})
ylabel('% change in Esr1+ spots/image')
yline(100, '--k')
yticks(0:50:200)
ylim([0 200])
axis square; set(gca, 'TickDir', 'out'); box off; grid off
pval = ranksum(con_norm, exp_norm);
effectsize = (mean(con_norm, 'omitnan') - mean(exp_norm, 'omitnan'))...
    ./std([con_norm; exp_norm]);
title(['d rank sum p=' num2str(pval) ', effect size =' num2str(effectsize)])

%--------------------------------------------------------------------------
% 6e. Example rat behavioral sensitivity to the blocks
%--------------------------------------------------------------------------
nexttile
exrat = 5; %
pre_hi_rat = ITIbyBlock_shEsr1_screen.('predox').high(exrat);
pre_hierr_rat = ITIbyBlock_shEsr1_screen.('predox').high_err(exrat);
pre_low_rat = ITIbyBlock_shEsr1_screen.('predox').low(exrat);
pre_lowerr_rat = ITIbyBlock_shEsr1_screen.('predox').low_err(exrat);
shEsr1_hi_rat = ITIbyBlock_shEsr1_screen.('duringdox').high(exrat);
shEsr1_hierr_rat = ITIbyBlock_shEsr1_screen.('duringdox').high_err(exrat);
shEsr1_low_rat = ITIbyBlock_shEsr1_screen.('duringdox').low(exrat);
shEsr1_lowerr_rat = ITIbyBlock_shEsr1_screen.('duringdox').low_err(exrat);
l = shadedErrorBar([1 2], [pre_low_rat pre_hi_rat],...
    [pre_lowerr_rat pre_hierr_rat], 'lineProps', {'-', 'Color',...
    doxstatecolors{1}, 'LineWidth', 0.5}); hold on
arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), l)
l = shadedErrorBar([1 2], [shEsr1_low_rat shEsr1_hi_rat],...
    [shEsr1_lowerr_rat shEsr1_hierr_rat], 'lineProps', {'-', 'Color',...
    doxstatecolors{2}, 'LineWidth', 0.5}); hold on
arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), l)
xlim([0.5 2.5])
xticks(1:2)
xticklabels({'Low', 'High'})
xlabel('Reward block')
yticks(-3:1:3)
ylim([-2.1 1.7])
ylabel('Detrended init.')
axis square; set(gca, 'TickDir', 'out'); box off
title('e Example rat')

%--------------------------------------------------------------------------
% 6f. Population level shEsr1 affect on behavioral sensitivity to the
% blocks
%--------------------------------------------------------------------------
nexttile
pre_delta = ITIbyBlock_shEsr1_screen.('predox').delta;
shEsr1_delta = ITIbyBlock_shEsr1_screen.('duringdox').delta;
pval_delta = signrank(pre_delta, shEsr1_delta);
effectsize_delta = effsize(pre_delta, shEsr1_delta);
%compare pre and post dox high and low
pre_high = ITIbyBlock_shEsr1_screen.('predox').high;
pre_low = ITIbyBlock_shEsr1_screen.('predox').low;
shEsr1_high = ITIbyBlock_shEsr1_screen.('duringdox').high;
shEsr1_low = ITIbyBlock_shEsr1_screen.('duringdox').low;
for rat = 1:length(ratlist_shEsr1)
    plot([1 2], [pre_low(rat) pre_high(rat)], '-',...
        linewidth=0.5, color = [0 0 0 0.2]); hold on
    plot([3 4], [shEsr1_low(rat) shEsr1_high(rat)], '-',...
        linewidth=0.5, color = [0 0 0 0.2]); hold on
end
plot(1, median(pre_low, 'omitnan'), '.k', markersize = 30); hold on
plot(2, median(pre_high, 'omitnan'), '.k', markersize = 30); hold on
plot(3, median(shEsr1_low, 'omitnan'), '.k', markersize = 30); hold on
plot(4, median(shEsr1_high, 'omitnan'), '.k', markersize = 30); hold on
errorbar(1, median(pre_low, 'omitnan'), sem(pre_low), 'k', linewidth=0.5, ...
    capsize=10); hold on
errorbar(2, median(pre_high, 'omitnan'), sem(pre_high), 'k', linewidth=0.5,...
    capsize=10); hold on
errorbar(3, median(shEsr1_low, 'omitnan'), sem(shEsr1_low), 'k', linewidth=0.5,...
    capsize=10); hold on
errorbar(4, median(shEsr1_high, 'omitnan'), sem(shEsr1_high), 'k', linewidth=0.5,...
    capsize=10); hold on
xlim([0.5 4.5])
xticks(1:1:4)
xticklabels({'Pre Low', 'Pre High', 'shEsr1 Low', 'shEsr1 High'})
yticks(-3:1:3)
ylim([-2.1 1.7])
ylabel('Detrended init.')
axis square; set(gca, 'TickDir', 'out'); box off
title(['f N=' num2str(length(ratlist_shEsr1))])
subtitle(['sign rank p = ' num2str(pval_delta) ', d = ' num2str(effectsize_delta)])

%--------------------------------------------------------------------------
% 6g. Model simulation of reduced gain on behavioral sensitivity to the
% reward blocks
%--------------------------------------------------------------------------
nexttile
for jj=1:2
    %manipulating RPE gain
    L = Ls{jj};
    Block = ratTrial.block;
    iti_block_avg = NaN(1,2);
    for bl = 2:3
        iti_block_avg(bl) = mean(L(Block==bl), 'omitnan');
    end    
    %low and high initiation time
    plot([1 2], [iti_block_avg(3) iti_block_avg(2)], '-',...
        Linewidth=0.5, color=doxstatecolors{jj}); hold on
    xlim([.5 2.5]);
    yticks([])
    ylabel('Modeled init. time');
    xticks(1:2)
    xticklabels({'Low'; 'High'});
    axis square; set(gca, 'TickDir', 'out'); box off
end
title('g')

%--------------------------------------------------------------------------
% 6h. Correlation of osmolality and estradiol in serum
%--------------------------------------------------------------------------
nexttile
plot(estradiol, RatConc, '.k', 'MarkerSize', 10); hold on
alpha(0.25)
xlabel('Estradiol (pg/ml)')
ylabel('Osmolality (mOsm/kgH2O)')
yticks(285:10:315)
bestfitline = lsline;
bestfitline.LineWidth = 0.5;
[R,p] = corrcoef(estradiol, RatConc, 'Rows', 'complete');
subtitle(['Pearson correlation: R=' num2str(R(1, 2)) ', p=' num2str(p(1, 2))])
axis square; grid off; set(gca, 'TickDir', 'out'); box off
title('h')

%--------------------------------------------------------------------------
% 6i. Population level shEsr1 affect on volume consumed
%--------------------------------------------------------------------------
%compare pre and post volume consumed
nexttile
pre_pro_volcon = ITIbyBlock_shEsr1_noscreen.('predox').('Proestrus').VolumeOverDuration;
pre_di_volcon = ITIbyBlock_shEsr1_noscreen.('predox').('Diestrus').VolumeOverDuration;
shEsr1_pro_volcon = ITIbyBlock_shEsr1_noscreen.('duringdox').('Proestrus').VolumeOverDuration;
shEsr1_di_volcon = ITIbyBlock_shEsr1_noscreen.('duringdox').('Diestrus').VolumeOverDuration;
for rat = 1:length(ratlist_shEsr1)
    plot([1 2], [pre_pro_volcon(rat) pre_di_volcon(rat)], '-',...
        linewidth=0.5, color = [0 0 0 0.2]); hold on
    plot([3 4], [shEsr1_pro_volcon(rat) shEsr1_di_volcon(rat)], '-',...
        linewidth=0.5, color = [0 0 0 0.2]); hold on    
end
plot(1, median(pre_pro_volcon, 'omitnan'), '.k', markersize = 30); hold on
plot(2, median(pre_di_volcon, 'omitnan'), '.k', markersize = 30); hold on
plot(3, median(shEsr1_pro_volcon, 'omitnan'), '.k', markersize = 30); hold on
plot(4, median(shEsr1_di_volcon, 'omitnan'), '.k', markersize = 30); hold on
errorbar(1, median(pre_pro_volcon, 'omitnan'), sem(pre_pro_volcon), 'k',...
    linewidth=0.5, capsize=10); hold on
errorbar(2, median(pre_di_volcon, 'omitnan'), sem(pre_di_volcon), 'k',...
    linewidth=0.5, capsize=10); hold on
errorbar(3, median(shEsr1_pro_volcon, 'omitnan'), sem(shEsr1_pro_volcon), 'k',...
    linewidth=0.5, capsize=10); hold on
errorbar(4, median(shEsr1_di_volcon, 'omitnan'), sem(shEsr1_di_volcon), 'k',...
    linewidth=0.5, capsize=10); hold on
xlim([0.5 4.5])
xticks(1:4)
xticklabels({'Pre Pro', 'Pre Di', 'shEsr1 Pro', 'shEsr1 Di'})
yticks(-1:0.2:3)
ylim([0.4 1.4])
ylabel(' Water consumed (μL)/session duration (s)')
axis square; set(gca, 'TickDir', 'out'); box off
pre_diff = pre_pro_volcon - pre_di_volcon;
shEsr1_diff = shEsr1_pro_volcon - shEsr1_di_volcon;
pval_volcon = signrank(pre_diff, shEsr1_diff);
effsize_volcon = effsize(pre_diff, shEsr1_diff);
title(['i N=' num2str(length(ratlist_shEsr1))])
subtitle(['sign rank p = ' num2str(pval_volcon)...
    ', effect size=' num2str(effsize_volcon)])


end