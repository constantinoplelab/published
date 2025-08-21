function PlotExtendedData6(datadir, codedir)
%PlotExtendedData6 - Plots Extended Data 4. 
% INPUTS:
%   datadir - Local directory where ProcessData_ExtendedData6.mat was saved after running ProcessData_ExtendedData6
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
load([datadir, 'ProcessData_ExtendedData6'],'Alignments', 'rewards', ...
    'da_rats_cpu','T','Stages','DA_G027','DA_err_G027', ...
    'high_cpu', 'low_cpu', 'err_high', 'err_low', 'AUC_byrat_CPu',...
    'stageffect_CPu', 'CPu_ratlist');
load([datadir, 'ProcessData_Figure2'], 'stageffect')

%% Create variables
mycolors = {'#3852A3', '#4A3F98', '#7B287C', '#BD1F43', '#EC2024'};
rew_vols = {'4', '8', '16', '32', '64'};

%% PLOT %%
figure;
set(gcf, 'Color', [1 1 1], 'Units', 'Inches', 'Position', [0, 0, 17, 9],...
    'PaperUnits', 'Inches', 'PaperSize', [9, 5], 'renderer','painters')

%--------------------------------------------------------------------------
% ED6d. Event aligned DA response, split by reward in mixed block for off-target rats.
%--------------------------------------------------------------------------
for a = 1:length(Alignments)
    DAplots = nan(1,length(rewards));
    nexttile
    for rew = 1:length(rewards)
        da_reward = NaN(length(CPu_ratlist), length(T));
        for rat = 1:length(CPu_ratlist)
            da_rat = da_rats_cpu{rat};
            da_reward(rat, :) = da_rat{rew, a};
        end
        %plot traces
        % subplot(2, length(Alignments), a)
        y = da_reward;
        DAplots(rew) = plot(T, mean(y, 'omitnan'),...
            'Color', mycolors{rew}); hold on
        sem = std(y, 'omitnan')./sqrt(size(y, 1));
        err = [mean(y, 'omitnan')-sem fliplr(mean(y, 'omitnan')+sem)];
        err(isnan(err)) = 0;
        h = fill([T fliplr(T)], err, 'k', 'FaceColor', mycolors{rew},...
            'LineStyle', 'none'); hold on
        set(h, 'facealpha', 0.25);
        set(gcf, 'Color', [1 1 1]);
        xlim([-0.5 1.25])
        set(gca, 'TickDir', 'out'); box off;
        xlabel(['Time from ', Alignments{a}]);
        axis square
    end
    ylim([-1 2.5])
    yticks(-2:1:2)
    line([0 0], ylim, 'Color', [0 0 0], 'LineStyle', '--');
    if a==length(Alignments)
        legend(DAplots, rew_vols, 'location', 'best')
    end
    title(['d ' Alignments{a}])
end

%--------------------------------------------------------------------------
% ED6e. Response to offer cue for all trials in mixed blocks, separated by reward volume and
% stage group for example rat, G027, baseline-corrected using the 0.05 to
% 0 s before offer cue.
%--------------------------------------------------------------------------
for e=1:length(Stages)
    nexttile
    for rew=1:length(rewards)
        da_mat = DA_G027{rew, e};
        da_err = DA_err_G027{rew, e};
        plot(T, da_mat, Color=mycolors{rew}, LineWidth=0.5,...
            displayname=rew_vols{rew}); hold on
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
    ylim([-1.2 2])
    xline(0, '--k');
    axis square
    title(['e Rat G027 ' Stages{e}])
end

%--------------------------------------------------------------------------
% ED6f. Response to offer cue, separated by reward block and
% stage group for example rat, G027. Gray box represents window used to
% calculate change in area under curve, 0 to 0.5 s delta AUC in g.
% Data is baseline-corrected using the 0.05 to 0 s before offer cue.
%--------------------------------------------------------------------------
cList = {'#FF0000', '#0000FF'}; %red, blue
for e=1:length(Stages)
    nexttile
    plot(T, high_cpu{e, 1}, color=cList{1},...
        linewidth=0.5, displayname='High'); hold on
    h = fill([T fliplr(T)], err_high{e}, 'k', 'FaceColor', cList{1},...
        'LineStyle', 'none');
    set(h, 'facealpha', 0.25);
    plot(T, low_cpu{e, 1}, color=cList{2},...
        linewidth=0.5, displayname='Low'); hold on
    h = fill([T fliplr(T)], err_low{e}, 'k', 'FaceColor', cList{2},...
        'LineStyle', 'none');
    set(h, 'facealpha', 0.25);
    ylabel('\Delta F/F')
    xlabel('Time from offer cue (s)')
    set(gca,'TickDir','out'); box off
    xlim([-0.5 1.25]);
    xticks(-1:0.5:1)
    axis square
    ylim([-0.8 1]); %set y-axis limit for all subplots
    yticks(-0.8:0.4:1.2)
    line([0 0], ylim, 'Color', [0 0 0], 'LineStyle', '--');
    xline([0.5 0.5], '--k'); hold on %end of AUC window
    title(['f Rat G027 ' Stages{e}]);
end

%--------------------------------------------------------------------------
% ED6g. Population average response to each reward volume, separated by reward
% block, min-max normalized, baseline-corrected using the 0.05 to 0 s
% before offer cue, *p<0.05.
%--------------------------------------------------------------------------
nexttile
linetypes = {'-','--'};
blocks = 1:3;
pvals = NaN(5, length(blocks));
effectsizes = NaN(5, length(blocks));
rewardsbyblock = {[1 2 3 4 5];[3 4 5];[1 2 3]};
mycolors = {'#7B287C', '#FF0000', '#0000FF'};
sem = @(xx) std(xx, 'omitnan') ./ sqrt(size(xx,1));
for bl=1:length(blocks)
    proAUC = NaN(length(CPu_ratlist),5);
    diAUC = NaN(length(CPu_ratlist),5);
    for e=1:length(Stages)
        clear thisplot        
        overrats = NaN(length(CPu_ratlist),5);
        for rat=1:length(CPu_ratlist)
            ratdata = AUC_byrat_CPu{rat};
            overrats(rat,:) = ratdata.(Stages{e}).AUC(bl, :);
        end
        if e==1
            proAUC = overrats;
        else
            diAUC = overrats;
        end
        thisplot = shadedErrorBar(rewards, median(overrats,'omitnan'),...
            sem(overrats), 'lineprops',{linetypes{e},...
            'color', mycolors{bl}, 'linewidth', 0.5}); hold on
        set(thisplot.edge,'LineStyle','none')
        arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), thisplot)
    end
    theserewards = rewardsbyblock{bl};
    for rew = theserewards(1):theserewards(end)
        pvals(rew,bl) = signrank(proAUC(:,rew),diAUC(:,rew));
        all_AUC = [proAUC(:,rew) diAUC(:,rew)];
        effectsizes(rew,bl) = (median(proAUC(:,rew)) -...
            median(diAUC(:,rew)))/std(all_AUC(:));
    end
end
title('h Proestrus -, Diestrus --')
axis square
set(gca, 'TickDir', 'out'); box off
xticks(rewards)
xlabel('Reward volume')
ylabel('AUC normalized')
yticks(-1:0.2:1)
disp('p-values:')
disp(pvals)
disp('effect sizes:')
disp(effectsizes)

%--------------------------------------------------------------------------
% ED6h. Histogram of stage effect (fertile-non-fertile) on deltaAUC for offer
% cue response to block (high - low)
%--------------------------------------------------------------------------
nexttile
stageffect_NAcc = stageffect; %loaded from Figure 2
histogram(stageffect_CPu, FaceColor='y',...
    binwidth=0.01); hold on
histogram(stageffect_NAcc, FaceColor='k',...
    binwidth=0.01); hold on
xline(0, '--k',LineWidth=0.5)
xline(median(stageffect_CPu), '-y',LineWidth=1)
xline(median(stageffect_NAcc), '-k',LineWidth=1)
xlabel('Proestrus - diestrus \Delta AUC')
ylabel('# Rats')
ylim([0 2.1]) %matches Fig 2
yticks(0:1:2)
set(gcf, 'Color', [1 1 1]);
set(gca, 'TickDir', 'out'); box off;
axis square
xlim([-0.16 0.16]) %matches Fig 2
pval_stageffect_CPu = signrank(stageffect_CPu);
title(['g N=' num2str(length(CPu_ratlist))...
    ' CPu sign rank diff from zero: p=' num2str(pval_stageffect_CPu)])

