function PlotExtendedData10(datadir, codedir)
%PlotExtendedData10 - Plots Extended Data 4. 
% INPUTS:
%   datadir - Local directory where ProcessData_ExtendedData10.mat was saved after running ProcessData_ExtendedData10
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
load([datadir, 'ProcessData_Figure1'], 'f_ratlist', 'ITIbyBlock')
load([datadir, 'ProcessData_ExtendedData10'], 'RatConc', 'RatGroups', ...
    'groups')


%% PLOT %%
figure;
set(gcf, 'Color', [1 1 1], 'Units', 'Inches', 'Position', [0, 0, 17, 9],...
    'PaperUnits', 'Inches', 'PaperSize', [9, 5], 'renderer','painters')

%--------------------------------------------------------------------------
%ED10a. Median intiation time in proestrus and diestrus stages
%--------------------------------------------------------------------------
nexttile
frats = length(f_ratlist);
cycle = {'Proestrus', 'Diestrus'};
for rat = 1:frats
    plot([1 2], [ITIbyBlock.(cycle{1}).ITIs_raw(rat)...
        ITIbyBlock.(cycle{2}).ITIs_raw(rat)],...
        '-', Color = [0 0 0 0.2], linewidth=0.5); hold on
end
plot([1 2], [median(ITIbyBlock.(cycle{1}).ITIs_raw, 'omitnan')...
    median(ITIbyBlock.(cycle{2}).ITIs_raw, 'omitnan')],...
    '.k', markersize=30); hold on
errorbar(1, median(ITIbyBlock.(cycle{1}).ITIs_raw, 'omitnan'),...
    sem(ITIbyBlock.(cycle{1}).ITIs_raw),...
    'k', capsize=10, linewidth=0.5); hold on
errorbar(2, median(ITIbyBlock.(cycle{2}).ITIs_raw, 'omitnan'),...
    sem(ITIbyBlock.(cycle{2}).ITIs_raw),...
    'k', capsize=10, linewidth=0.5); hold on
xlim([0.5 2.5])
xticks([1 2])
xticklabels(cycle)
ylabel('Initiation time (s)')
ylim([0.5 8.2])
yticks(0:2:8)
axis square; grid off; set(gca, 'TickDir', 'out'); box off
pval = signrank(ITIbyBlock.(cycle{1}).ITIs_raw,...
    ITIbyBlock.(cycle{2}).ITIs_raw);
effectsize = effsize(ITIbyBlock.(cycle{1}).ITIs_raw,...
    ITIbyBlock.(cycle{2}).ITIs_raw);
subtitle(['sign rank p=' num2str(pval) ', effect size = ' num2str(effectsize)])
title('a')

%--------------------------------------------------------------------------
%ED10b. Number of trials (corrected for session duration) in proestrus and 
% diestrus stages
%--------------------------------------------------------------------------
nexttile
for rat = 1:frats
    plot([1 2], [ITIbyBlock.(cycle{1}).TrialsOverDuration(rat)...
        ITIbyBlock.(cycle{2}).TrialsOverDuration(rat)],...
        '-', Color = [0 0 0 0.2], linewidth=0.5); hold on
end
plot([1 2], [median(ITIbyBlock.(cycle{1}).TrialsOverDuration, 'omitnan')...
    median(ITIbyBlock.(cycle{2}).TrialsOverDuration, 'omitnan')],...
    '.k', markersize=30); hold on
errorbar(1, median(ITIbyBlock.(cycle{1}).TrialsOverDuration, 'omitnan'),...
    sem(ITIbyBlock.(cycle{1}).TrialsOverDuration),...
    'k', capsize=10, linewidth=0.5); hold on
errorbar(2, median(ITIbyBlock.(cycle{2}).TrialsOverDuration, 'omitnan'),...
    sem(ITIbyBlock.(cycle{2}).TrialsOverDuration),...
    'k', capsize=10, linewidth=0.5); hold on
xlim([0.5 2.5])
xticks([1 2])
xticklabels(cycle)
ylabel('Number of trials/session duration')
yticks(0:0.1:0.3)
ylim([0.02 0.27])
grid off; axis square; set(gca, 'TickDir', 'out'); box off
pval = signrank(ITIbyBlock.(cycle{1}).TrialsOverDuration,...
    ITIbyBlock.(cycle{2}).TrialsOverDuration);
effectsize = effsize(ITIbyBlock.(cycle{1}).TrialsOverDuration,...
    ITIbyBlock.(cycle{2}).TrialsOverDuration);
subtitle(['sign rank p=' num2str(pval) ', effect size = ' num2str(effectsize)])
title('b')

%--------------------------------------------------------------------------
%ED10c. Volume of water consumed (corrected for session duration) in 
% proestrus and diestrus stages
%--------------------------------------------------------------------------
nexttile
for rat = 1:frats
    plot([1 2], [ITIbyBlock.(cycle{1}).VolumeOverDuration(rat)...
        ITIbyBlock.(cycle{2}).VolumeOverDuration(rat)],...
        '-', Color = [0 0 0 0.2], linewidth=0.5); hold on
end
plot([1 2], [median(ITIbyBlock.(cycle{1}).VolumeOverDuration, 'omitnan')...
    median(ITIbyBlock.(cycle{2}).VolumeOverDuration, 'omitnan')],...
    '.k', markersize=30); hold on
errorbar(1, median(ITIbyBlock.(cycle{1}).VolumeOverDuration, 'omitnan'),...
    sem(ITIbyBlock.(cycle{1}).VolumeOverDuration),...
    'k', capsize=10, linewidth=0.5); hold on
errorbar(2, median(ITIbyBlock.(cycle{2}).VolumeOverDuration, 'omitnan'),...
    sem(ITIbyBlock.(cycle{2}).VolumeOverDuration),...
    'k', capsize=10, linewidth=0.5); hold on
xlim([0.5 2.5])
xticks([1 2])
xticklabels(cycle)
ylabel('Volume consumed/session duration')
ylim([0.3 2.5])
yticks(0.4:0.4:3.2)
grid off; axis square; set(gca, 'TickDir', 'out'); box off
pval = signrank(ITIbyBlock.(cycle{1}).VolumeOverDuration,...
    ITIbyBlock.(cycle{2}).VolumeOverDuration);
effectsize = effsize(ITIbyBlock.(cycle{1}).VolumeOverDuration,...
    ITIbyBlock.(cycle{2}).VolumeOverDuration);
subtitle(['sign rank p=' num2str(pval) ', effect size = ' num2str(effectsize)])
title('c')

%--------------------------------------------------------------------------
%ED10d. Median serum osmolality by group
%--------------------------------------------------------------------------
nexttile
outdata = cell(length(groups), 1);
for g = 1:length(groups) %edit: remove males
    thisgroup = groups(g);
    group_conc = RatConc(strcmp(RatGroups, thisgroup));
    y = median(group_conc, 'omitnan');
    this_sem = std(group_conc, 'omitnan')./sqrt(sum(~isnan(group_conc)));
    plot(g, y, '.k', markersize=30); hold on
    errorbar(g, y, this_sem, 'k', 'LineWidth', 0.5, capsize=10)
    outdata{g} = group_conc;
end
xticks(1:length(groups))
xticklabels(groups)
ylabel('Osmolality (mOsm/kgH2O)')
yticks(200:5:400)
xlim([0.5 5.5])
axis square; grid off; set(gca, 'TickDir', 'out'); box off
pval = kruskalwallis([outdata{1}; outdata{2}; outdata{3}; outdata{4};...
    outdata{5}], [ones(length(outdata{1}), 1);...
    ones(length(outdata{2}), 1)*2;...
    ones(length(outdata{3}), 1)*3;...
    ones(length(outdata{4}), 1)*4;...
    ones(length(outdata{5}), 1)*5], 'off');
effectsize = effsize_ind(outdata{2},outdata{4});
title(['d kruskal wallis p=' num2str(pval)])
subtitle(['late p vs d effect size (independent) = ' num2str(effectsize)])

end



