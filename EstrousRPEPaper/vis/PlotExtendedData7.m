function PlotExtendedData7(datadir, codedir)
%PlotExtendedData7 - Plots Extended Data 4. 
% INPUTS:
%   datadir - Local directory where ProcessData_ExtendedData7.mat was saved after running ProcessData_ExtendedData7
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
load([datadir, 'ProcessData_ExtendedData7'],'NAcc_ratlist','numbins', ...
    'da_rats_SideOn','T', 'bins', 'delay_legend', 'da_rats_SideOff', ...
    'Stages', 'da_rats_stages', 'AUC_stages')

%% Create variables
cycle_colors = {'#E87003';'#7358C6'}; %dark orange, purple

%% PLOT %%
figure;
set(gcf, 'Color', [1 1 1], 'Units', 'Inches', 'Position', [0, 0, 17, 9],...
    'PaperUnits', 'Inches', 'PaperSize', [9, 5], 'renderer','painters')

%--------------------------------------------------------------------------
% ED7a. Dopamine during the delay, split by delay length.
% Includes negative ramping, which suggests state inference.
%--------------------------------------------------------------------------
subplot(1, 5, 1)
mycolors = gray(numbins+1);
mycolors = flipud(mycolors);
da_bins_rats = cell(length(bins), 1);
for j = 2:length(bins)
    da_bin_rats = NaN(length(NAcc_ratlist), size(T, 2));
    for rat = 1:length(NAcc_ratlist)
        darat = da_rats_SideOn{rat};
        da_bin_rats(rat, :) = darat{j};
    end
    da_bins_rats{j} = da_bin_rats;
end
DAplots = nan(1,length(bins));
for j = 2:length(bins)
    DA = da_bins_rats{j};
    y = DA;
    DAplots(j) = plot(T, mean(y, 'omitnan'),...
        'Color', mycolors(j,:), 'DisplayName', delay_legend{j}); hold on
    thiserr = std(y, 'omitnan')./sqrt(size(y, 1));
    thisplot = shadedErrorBar(T, mean(y, 'omitnan'), thiserr, 'lineprops',{'-',...
        'color', mycolors(j, :), 'linewidth', 0.5});
    arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), thisplot)
end
hleg = legend(DAplots(2:end), delay_legend(2:end), 'location', 'best');
title(hleg,' delay (less than or equal to, s)')
xlim([-0.75 6]);
xticks(1:1:6)
xlabel('Time since delay start (s)');
axis square; set(gca, 'TickDir', 'out'); box off
xline(0, '--k', HandleVisibility='off')
title('a')

%--------------------------------------------------------------------------
% ED7b. Dopamine at the time of the side light turning off (cues that the
% delay is over and the reward is available), split by delay length.NAcc_ratlist
%--------------------------------------------------------------------------
subplot(1, 5, 2)
da_bins_rats = cell(length(bins), 1);
for j = 2:length(bins)
    da_bin_rats = NaN(length(NAcc_ratlist), size(T, 2));
    for rat = 1:length(NAcc_ratlist)
        darat = da_rats_SideOff{rat};
        da_bin_rats(rat, :) = darat{j};
    end
    da_bins_rats{j} = da_bin_rats;
end
for j = 2:length(bins)
    DA = da_bins_rats{j};
    y = DA;
    plot(T, mean(y, 'omitnan'),...
        'Color', mycolors(j,:), 'DisplayName', delay_legend{j}); hold on
    thiserr = std(y, 'omitnan')./sqrt(size(y, 1));
    thisplot = shadedErrorBar(T, mean(y, 'omitnan'), thiserr, 'lineprops',{'-',...
        'color', mycolors(j, :), 'linewidth', 0.5});
    arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), thisplot)
end
set(gcf, 'Color', [1 1 1]);
xlim([-0.75 1.5]);
xlabel('Time since delay end (s)');
ylim([-0.1 1.6])
axis square; set(gca, 'TickDir', 'out'); box off
xline(0, '--k', HandleVisibility='off')
title('b')

%--------------------------------------------------------------------------
% ED7c. Phasic response at the time of the side light turning off, split by 
% estrous stage (additional evidence that RPE is enhanced in proestrus).
% Time course and quantification of area under the curve.
%--------------------------------------------------------------------------
for s = 1:length(Stages)
    da_bins_rats = cell(length(bins), 1);
    for j = 2:length(bins)
        da_bin_rats = NaN(length(NAcc_ratlist), size(T, 2));
        for rat = 1:length(NAcc_ratlist)
            da_rat = da_rats_stages{rat};
            da_bin_rats(rat, :) = da_rat{s, j};
        end
        da_bins_rats{j} = da_bin_rats;
    end
    for j = 2:length(bins)
        subplot(1, 5, 2+s);
        DA = da_bins_rats{j};
        y = DA;
        thisplot = shadedErrorBar(T, mean(y, 'omitnan'), sem(y), ...
            'lineprops',{'-','color', mycolors(j, :), 'linewidth', 0.5});
        arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), ...
            line.edge), thisplot)
    end
    xlim([-0.5 1]);
    xlabel('Time since delay end (s)');    
    axis square; set(gca, 'TickDir', 'out'); box off
    title(['c ' Stages{s}])
end
% set y-axis limit for all subplots
for s = 1:length(Stages)
    subplot(1, 5, 2+s)
    ylim([-0.2 1.65]);
    yticks(0:1:2)
    xline(0, '--k', HandleVisibility='off')
end

%--------------------------------------------------------------------------
% ED7d. AUCs of average
%--------------------------------------------------------------------------
subplot(1,5,5)
for s=1:length(Stages)
    AUC_thisstage = NaN(length(bins), length(NAcc_ratlist));
    for rat = 1:length(NAcc_ratlist)
        AUC_rat = AUC_stages{rat};
        AUC_thisstage(:, rat) = AUC_rat(s, :)';
    end
    AUC_thisstage = AUC_thisstage';
    thisplot = shadedErrorBar(bins(1:end-1), median(AUC_thisstage(:, 2:end), 'omitnan'),...
        sem(AUC_thisstage(:, 2:end)), 'lineprops',{'-',...
            'color', cycle_colors{s}, 'linewidth', 0.5}); hold on; drawnow
    arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), thisplot)
    if s==1
        AUC_pro = AUC_thisstage;
    elseif s==2
        AUC_di = AUC_thisstage;
    end
end
xlim([bins(1)-0.3 bins(end-1)+0.3])
xlabel('Reward delay (s)');  
ylabel('Delay end AUC');
ylim([0.1 0.7])
yticks(0.2:0.2:0.8)
xticks(0:1:8)
axis square; set(gca, 'TickDir', 'out'); box off
%%Stats
pvals = length(bins)-1;
for bin=2:length(bins)
    pvals(bin-1) = signrank(AUC_pro(:, bin), AUC_di(:, bin));
end
disp(pvals)
title('d')

end



