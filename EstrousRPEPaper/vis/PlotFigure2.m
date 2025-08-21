function PlotFigure2(datadir, codedir)
%PlotFigure2 - Plots Figure 2. 
% INPUTS:
%   datadir - Local directory where ProcessData_Figure2.mat was saved after running ProcessData_Figure2
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
load([datadir, 'ProcessData_Figure2'],...
    'NAcc_ratlist', 'Stages', 'times_resample',...
    'sorted_data_rew_rat', 'reward_da_rats',...
    'T', 'sorted_data_16', 'high_rats', 'low_rats', 'high_stages_rat',...
    'low_stages_rat', 'err_high_stages_rat', 'err_low_stages_rat',...
    'stageffect', 'reward_da_stages_rat', 'reward_da_err_stages_rat',...
    'AUC_byrat', 'delta_good', 'E2_good', 'stages');

%% PLOT %%
figure;
set(gcf, 'Color', [1 1 1], 'Units', 'Inches', 'Position', [0, 0, 17, 9],...
    'PaperUnits', 'Inches', 'PaperSize', [9, 5], 'renderer','painters')

%--------------------------------------------------------------------------
%2c.  Heatmap of motion-corrected GRABDA signal in response to
% reward offer cue across over 16,000 trials recorded from example rat,
% G024, sorted by reward volume.
%--------------------------------------------------------------------------
nexttile
ratname = 'G024';
imagesc(times_resample, 1:size(sorted_data_rew_rat,1),...
    sorted_data_rew_rat(:,1:end-2)); hold on
colorbar
set(gca, 'TickDir', 'out'); box off
xlim([min(times_resample) max(times_resample)]);
clim([-2 2]);
xlabel('Time from offer cue (s)');
ylabel('Trial');
xlim([-0.5 1.25])
xticks(-1:0.5:1)
xline([0 0], '-k');
yticks([1 16000])
xline([0 0], '-k', linewidth=0.5);
title(['c Rat ' ratname])

%--------------------------------------------------------------------------
%2d. Average response to offer cue across the 9 rats recorded from,
% separated by reward volume.
%--------------------------------------------------------------------------
nexttile
rewards = [1 2 3 4 5];
rewardnames = {'4', '8', '16', '32', '64'};
mycolors = {'#3852A3', '#4A3F98', '#7B287C', '#BD1F43', '#EC2024'};
for rew = 1:length(rewards)
    da_reward = NaN(length(NAcc_ratlist), length(T));
    for rat = 1:length(NAcc_ratlist)
        da_rat = reward_da_rats{rat};
        da_reward(rat, :) = da_rat{rew, 1};
    end
    y = da_reward;
    plot(T, mean(y, 'omitnan'), color=mycolors{rew},...
        DisplayName=rewardnames{rew}, linewidth=0.5); hold on
    sem_rew = std(y, 'omitnan')./sqrt(size(y, 1));
    err = [mean(y, 'omitnan')-sem_rew fliplr(mean(y, 'omitnan')+sem_rew)];
    err(isnan(err)) = 0;
    h = fill([T fliplr(T)], err, 'k', FaceColor=mycolors{rew},...
        LineStyle='none'); hold on
    set(h, 'facealpha', 0.25);
end
set(gcf, 'Color', [1 1 1]);
xlim([-0.5 1.25])
xticks(-1:0.5:1)
ylim([-0.8 2.1])
yticks(-2:1:2)
set(gca, 'TickDir', 'out'); box off;
xlabel('Time from offer cue (s)');
ylabel('\Delta F/F');
axis square
xline(0, '--k', 'linewidth', 0.5);
title(['d Mean over rats n = ', num2str(length(NAcc_ratlist))]);

%--------------------------------------------------------------------------
%2e.  Heatmap of motion-corrected GRABDA signal in response to reward
% offer cue for 16ul trials across over 7,000 trials recorded from example rat,
% G022, sorted by reward block.
%--------------------------------------------------------------------------
nexttile
ratname = 'G022';
imagesc(times_resample, 1:size(sorted_data_16,1), sorted_data_16(:,1:end-2)); hold on
yline(find(diff(sorted_data_16(:,end-1))>0), '--w', 'LineWidth', 0.5)
colorbar
set(gca, 'TickDir', 'out'); box off
xlim([min(times_resample) max(times_resample)]);
clim([-0.75 0.75]);
xlabel('Time from offer cue (s)');
ylabel('Trial');
xlim([-0.5 1.25])
xticks(-1:0.5:1)
yticks([1 3000])
xline([0 0], '-k', linewidth=0.5);
title(['e Rat ' ratname])

%--------------------------------------------------------------------------
%2f.  Average response to 16ul offer cue across the 9 rats recorded
% from, separated by reward block.
%--------------------------------------------------------------------------
nexttile
hi_event = NaN(length(NAcc_ratlist), size(T, 2));
lo_event = NaN(length(NAcc_ratlist), size(T, 2));
for rat = 1:length(NAcc_ratlist)
    hi = high_rats{rat};
    lo = low_rats{rat};
    hi_event(rat, :) = hi;
    lo_event(rat, :) = lo;
end
mean_hi_event = mean(hi_event, 'omitnan');
mean_lo_event = mean(lo_event, 'omitnan');
plot(T, mean_hi_event, 'r', 'LineWidth', 0.5, ...
    'DisplayName', 'High'); hold on
plot(T, mean_lo_event, 'b', 'LineWidth', 0.5,...
    'DisplayName', 'Low'); hold on
sem_hi = std(hi_event, 'omitnan')./sqrt(size(hi_event, 1));
err_hi = [mean_hi_event-sem_hi fliplr(mean_hi_event+sem_hi)];
err_hi(isnan(err_hi)) = 0;
h_hi = fill([T fliplr(T)], err_hi, 'r', 'LineStyle', 'none');
set(h_hi, 'facealpha', 0.25);
sem_lo = std(lo_event, 'omitnan')./sqrt(size(lo_event, 1));
err_lo = [mean_lo_event-sem_lo fliplr(mean_lo_event+sem_lo)];
err_lo(isnan(err_lo)) = 0;
h_lo = fill([T fliplr(T)], err_lo, 'b', 'LineStyle', 'none');
set(h_lo, 'facealpha', 0.25);
xline(0, '--k');
xlim([-0.5 1.25]);
xticks(-1:0.5:1)
yticks(-1:1:1)
ylim([-1.3 1])
xlabel('Time from offer cue (s)');
set(gcf, 'Color', [1 1 1]);
set(gca, 'TickDir', 'out'); box off; axis square
title(['f Mean over rats n = ', num2str(length(NAcc_ratlist))]);

%--------------------------------------------------------------------------
%2g. Response to offer cue, separated by reward block and
% stage group for example rat. Gray box represents window used to
% calculate change in area under curve, 0 to 0.5 s delta AUC in h.
% Data is baseline-corrected using the 0.05 to 0 s before offer cue.
%--------------------------------------------------------------------------
ratname = 'G016';
cList = {'#FF0000', '#0000FF'}; %red, blue
range = [0 0];
ax = NaN(length(Stages), 1);
for s=1:length(Stages)
    ax(s) = nexttile;
    plot(T, high_stages_rat{s, 1}, color=cList{1},...
        linewidth=0.5, displayname='High'); hold on
    h = fill([T fliplr(T)], err_high_stages_rat{s}, 'k', 'FaceColor', cList{1},...
        'LineStyle', 'none');
    set(h, 'facealpha', 0.25);
    plot(T, low_stages_rat{s, 1}, color=cList{2},...
        linewidth=0.5, displayname='Low'); hold on
    h = fill([T fliplr(T)], err_low_stages_rat{s}, 'k', 'FaceColor', cList{2},...
        'LineStyle', 'none');
    set(h, 'facealpha', 0.25);
    ylabel('\Delta F/F')
    xlabel('Time from offer cue (s)')
    set(gca,'TickDir','out'); box off
    xlim([-0.5 1.25]);
    xticks(-1:0.5:1)
    axis square
    yl = ylim;
    if yl(1) < range(1)
        range(1) = yl(1);
    end
    if yl(2) > range(2)
        range(2) = yl(2);
    end
    title(['g Rat ' ratname ': ' Stages{s}])
    ylim(ax(s),[-0.5 0.7])
    yticks(-0.8:0.4:0.8)
    xline(0.5, '--k'); hold on %end of AUC window
    xline(0, '--k');
end

%--------------------------------------------------------------------------
%2h. Histogram of stage effect (proestrus-diestrus) on deltaAUC for offer
% cue response, *\p<0.05, Wilcoxon signed-rank test for difference from zero,
%baseline corrected
%--------------------------------------------------------------------------
nexttile
histogram(stageffect, FaceColor='k',...
    binwidth=0.01); hold on
xline(median(stageffect), '--k',LineWidth=0.5)
xlabel('Proestrus - diestrus \Delta AUC')
ylabel('# rats')
yticks(0:1:2)
xline(0, '--k', 'LineWidth', 0.5)
set(gcf, 'Color', [1 1 1]);
set(gca, 'TickDir', 'out'); box off; axis square
xlim([-0.16 0.16])
ylim([0 2.1])
pval_stageeffect = signrank(stageffect);
title(['h p=' num2str(pval_stageeffect)])

%--------------------------------------------------------------------------
%2i. Response to offer cue as a function of estradiol
%--------------------------------------------------------------------------
nexttile
[R,P] = corrcoef(E2_good, delta_good);
plot(E2_good, delta_good, '.k'); hold on
lsline
plts = NaN(4,1);
plts(1) = plot(E2_good(strcmp(stages, 'Proestrus')),...
    delta_good(strcmp(stages, 'Proestrus')), '.',...
    Color='#E87003', MarkerSize=20); hold on
plts(2) = plot(E2_good(strcmp(stages, 'Estrus')),...
    delta_good(strcmp(stages, 'Estrus')), '.',...
    Color='#DD965B', MarkerSize=20); hold on
plts(3) = plot(E2_good(strcmp(stages, 'Metestrus')),...
    delta_good(strcmp(stages, 'Metestrus')), '.',...
    Color='#A495E5', MarkerSize=20); hold on
plts(4) = plot(E2_good(strcmp(stages, 'Diestrus')),...
    delta_good(strcmp(stages, 'Diestrus')), '.',...
    Color='#7358C6', MarkerSize=20); hold on
legend(plts, {'Proestrus', 'Estrus', 'Metestrus', 'Diestrus'})
ylim([0 0.8])
yticks(0:0.4:0.8)
xlabel('Estradiol (pg/ml)')
ylabel('\Delta AUC') %
title(['i R=' num2str(R(2,1)) ', p=' num2str(P(2,1))])
grid off; axis square; set(gca, 'TickDir', 'out'); box off

% xbin = linspace(min(E2_good), max(E2_good), 9);
% ybin = nan(1, length(xbin));
% ybiner = ybin;
% for j = 2:length(ybin)
%     these = find(E2_good>xbin(j-1) & E2_good<=xbin(j));
%     ybin(j) = mean(delta_good(these));
%     ybiner(j) = std(delta_good(these))./sqrt(length(these));
% end
% 
% centeredx = xbin-((xbin(2)-xbin(1))/2); %x-values centered for each bin
% subplot(1,2,2);
% l = shadedErrorBar(centeredx, ybin, ybiner, 'lineprops', '-k'); 
% arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), l)
% xlabel('Estradiol (pg/ml)');
% ylabel('\Delta AUC (high - low)'); 
% title(strcat([num2str(sum(~isnan(delta_good))) ' sessions']));
% grid off; axis square; set(gca, 'TIckDir', 'out'); box off;
% %ylim([.5 2.5]);
% xlim([30 200])
% set(gcf, 'Color', [1 1 1]);

%--------------------------------------------------------------------------
%2j. Response to offer cue for all trials in mixed blocks, separated by reward volume and
% stage group for example rat, G037, baseline-corrected using the 0.05 to
% 0 s before offer cue.
%--------------------------------------------------------------------------
ratname = 'G016';
range = [0 0];
ax = NaN(length(Stages), 1);
for s=1:length(Stages)
    ax(s) = nexttile;
    for rew=1:length(rewards)
        da_mat = reward_da_stages_rat{rew, s};
        da_err = reward_da_err_stages_rat{rew, s};
        plot(T, da_mat, Color=mycolors{rew}, LineWidth=0.5,...
            displayname=rewardnames{rew}); hold on
        h = fill([T fliplr(T)], da_err, 'k', facecolor=mycolors{rew},...
            LineStyle='none');
        set(h, 'facealpha', 0.25);
    end
    xlim([-0.5 1.25])
    xticks(-1:0.5:1)
    yticks(-1:1:2)
    set(gca, 'TickDir', 'out'); box off
    xlabel('Time from offer cue (s)');
    ylabel('\Delta F/F')
    yl = ylim;
    if yl(1) < range(1)
        range(1) = yl(1);
    end
    if yl(2) > range(2)
        range(2) = yl(2);
    end
    axis square
    title(['j Rat ' ratname ': ' Stages{s}]);
    xline(0, '--k');
    xline(0.5, '--k');
end
for s = 1:length(Stages)
    ylim(ax(s),range)
end

%--------------------------------------------------------------------------
%2k left. Population average response to each reward volume, separated by reward
% block, min-max normalized, baseline-corrected using the 0.05 to 0 s
% before offer cue, *p<0.05.
%--------------------------------------------------------------------------
nexttile
linetypes = {'-','--'};
blocks = 1:3;
allrewards = [4 8 16 32 64];
pvals = NaN(5, length(blocks));
effectsizes = NaN(5, length(blocks));
rewards = {[1 2 3 4 5];[3 4 5];[1 2 3]};
mycolors = {'#7B287C', '#FF0000', '#0000FF'};
sem = @(xx) std(xx, 'omitnan') ./ sqrt(sum(~isnan(xx)));
for bl=1:length(blocks)
    proAUC = NaN(length(NAcc_ratlist),5);
    diAUC = NaN(length(NAcc_ratlist),5);
    for s=1:length(Stages)
        clear thisplot        
        overrats = NaN(length(NAcc_ratlist),5);
        for rat=1:length(NAcc_ratlist)
            ratdata = AUC_byrat{rat};
            overrats(rat,:) = ratdata.(Stages{s}).AUC_norm(bl, :);
        end
        if s==1
            proAUC = overrats;
        else
            diAUC = overrats;
        end
        thisplot = shadedErrorBar(allrewards, median(overrats,'omitnan'),...
            sem(overrats), 'lineprops',{linetypes{s},...
            'color', mycolors{bl}, 'linewidth', 0.5}); hold on
        set(thisplot.edge,'LineStyle','none')
        arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), thisplot)
    end
    theserewards = rewards{bl};
    for rew = theserewards(1):theserewards(end)
        pvals(rew,bl) = signrank(proAUC(:,rew),diAUC(:,rew));
        all_AUC = [proAUC(:,rew) diAUC(:,rew)];
        effectsizes(rew,bl) = (median(proAUC(:,rew)) -...
            median(diAUC(:,rew)))/std(all_AUC(:));
    end
end
title('k Proestrus -, Diestrus --')
set(gcf, 'Color', [1 1 1]);
set(gca, 'TickDir', 'out'); box off; axis square
xticks(allrewards)
xlabel('Reward volume')
ylabel('AUC normalized')
ylim([-0.03 1.03])
yticks(0:0.2:1)
disp('p-values:')
disp(pvals)
disp('effect sizes:')
disp(effectsizes)

%--------------------------------------------------------------------------
%2k right top. Individual responses to high block largest reward, not normalized
%--------------------------------------------------------------------------
nexttile
overrats_64ul_high = NaN(length(Stages), length(NAcc_ratlist));
overrats_64ul_mixed = NaN(length(Stages), length(NAcc_ratlist));
for s=1:length(Stages)
    for rat=1:length(NAcc_ratlist)
        ratdata = AUC_byrat{rat};
        high = ratdata.(Stages{s}).AUC(2, :);
        mixed = ratdata.(Stages{s}).AUC(1, :);
        overrats_64ul_high(s, rat) = high(end);
        overrats_64ul_mixed(s, rat) = mixed(end);
    end
end
diff_high = overrats_64ul_high(1, :)-overrats_64ul_high(2, :);
these_binedges = [-0.4125:0.025:-0.0125 0.0125:0.025:0.4125];
histogram(diff_high, facecolor='k', binwidth=0.1, binedges=these_binedges)
xline(0, '--k')
xline(median(diff_high), '--k')
xlim([-0.4 0.4])
ylim([0 3.25])
yticks(0:1:3)
ylabel('# rats')
set(gcf, 'Color', [1 1 1]);
set(gca, 'TickDir', 'out'); box off; axis square
title('k High')

%--------------------------------------------------------------------------
%2k right bottom. Individual responses to mixed block largest reward, not normalized
%--------------------------------------------------------------------------
nexttile
diff_mix = overrats_64ul_mixed(1, :)-overrats_64ul_mixed(2, :);
these_binedges = [-0.4125:0.025:-0.0125 0.0125:0.025:0.4125];
histogram(diff_mix, facecolor='k', binwidth=0.1, binedges=these_binedges)
xline(0, '--k')
xline(median(diff_high), '--k')
xlim([-0.4 0.4])
ylim([0 3.25])
yticks(0:1:3)
ylabel('# rats')
set(gcf, 'Color', [1 1 1]);
set(gca, 'TickDir', 'out'); box off; axis square
title('k Mixed')

end