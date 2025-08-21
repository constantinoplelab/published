function PlotExtendedData2(datadir, codedir)
%PlotExtendedData1 - Plots Extended Data 2. 
% INPUTS:
%   datadir - Local directory where ProcessData_ExtendedData2.mat was saved after running ProcessData_ExtendedData2
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
load([datadir 'ProcessData_ExtendedData2'], 'WTbyBlock', 'f_ratlist');
cycle = {'Proestrus', 'Estrus', 'Metestrus', 'Diestrus'};

%% PLOT %%
figure;
set(gcf, 'Color', [1 1 1], 'Units', 'Inches', 'Position', [0, 0, 17, 9],...
    'PaperUnits', 'Inches', 'PaperSize', [9, 5], 'renderer','painters')

%--------------------------------------------------------------------------
%ED2a. Mean + SEM wait time by block and by stage for example rat
%--------------------------------------------------------------------------
nexttile
x = [1 2 3 4 5];
b_colors = [{[1 0 0]};{[0 0 1]}];
linetypes = {'-', '--'};
frats = length(f_ratlist);
for e = 1:length(cycle)

    hi = NaN(length(frats), length(x));
    lo = NaN(length(frats), length(x));

    for rat = 1: frats     
            hi(rat, :) = WTbyBlock.(cycle{e}).mean_high{rat}; %mean high wts per rat
            lo(rat, :) = WTbyBlock.(cycle{e}).mean_low{rat}; %mean low wts per rat
    end

    thisplot = shadedErrorBar(x, mean(hi, 'omitnan'), sem(hi), 'lineprops', {linetypes{e},....
        'Color', b_colors{1}, 'LineWidth', 0.5}); hold on
    set(thisplot.edge,'LineStyle','none')
    arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), thisplot)
    thisplot = shadedErrorBar(x, mean(lo, 'omitnan'), sem(lo), 'lineprops', {linetypes{e},....
        'Color', b_colors{2}, 'LineWidth', 0.5}); hold on
    set(thisplot.edge,'LineStyle','none')
    arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), thisplot)

    %find number of rats included in figure
    N = sum(~isnan(lo(:, 3)));

    xticks([1 2 3 4 5])        
    xlim([0.5 5.5])
    ylabel('Wait time (s)')
    xlabel('Reward volume')
    grid off
    set(gca, 'TickDir', 'out'); box off
    axis square    
end 
ylim([10.5 13.5])
yticks(10:1:14)
xticklabels([4 8 16 32 64])
title(['a N=', num2str(N)])

%--------------------------------------------------------------------------
%ED2b. Delta summary statistics
%--------------------------------------------------------------------------
nexttile

delta_pro = WTbyBlock.(cycle{1}).delta;
delta_est = WTbyBlock.(cycle{2}).delta;
delta_met = WTbyBlock.(cycle{3}).delta;
delta_di = WTbyBlock.(cycle{4}).delta;

%plot average and sem 
for e = 1:length(cycle)
    delta = WTbyBlock.(cycle{e}).delta;
    plot(e, median(delta, 'omitnan'), '.k', markersize=30,...
        LineWidth=0.5); hold on
    errorbar(e, median(delta, 'omitnan'), sem(delta), LineWidth=0.5,...
        color='k',capsize=10); hold on        
end
kw_pval = kruskalwallis([delta_pro'; delta_est'; delta_met'; delta_di'],...
    [ones(frats, 1); ones(frats, 1)*2; ones(frats, 1)*3; ones(frats, 1)*4],...
    'off');
delta_p = signrank(delta_pro, delta_di);
effectsize = effsize(delta_pro, delta_di);
xlim([0.5 length(cycle)+0.5]);
xticks(1:length(cycle))
xticklabels(cycle)
ylabel('Wait time for 16ul low - high')
ylim([1 1.9])
grid off
set(gca, 'TickDir', 'out'); box off
axis square
title(['b kruskal wallis p=' num2str(kw_pval)])
subtitle(['p vs d sign rank p=' num2str(delta_p) ', effect size = ' num2str(effectsize)])


end