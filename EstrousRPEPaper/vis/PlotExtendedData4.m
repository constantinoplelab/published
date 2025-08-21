function PlotExtendedData4(datadir, codedir)
%PlotExtendedData4 - Plots Extended Data 4. 
% INPUTS:
%   datadir - Local directory where ProcessData_ExtendedData4.mat was saved after running ProcessData_ExtendedData4
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
load([datadir 'ProcessData_ExtendedData4'], 'NAcc_ratlist', 'rewards', ...
    'Alignments', 'DAbyrat', 'DAerrbyrat', 'T')

%% Create variables
mycolors = {'#3852A3', '#4A3F98', '#7B287C', '#BD1F43', '#EC2024'}; %ROYGBV
rew_vols = {'4', '8', '16', '32', '64'};


%% PLOT %%
%--------------------------------------------------------------------------
% Event-aligned DA by reward for each rat. New figure per rat.
%--------------------------------------------------------------------------
for rat = 1:length(NAcc_ratlist)

    figure;
    set(gcf, 'Color', [1 1 1], 'Units', 'Inches', 'Position', [0, 0, 17, 9],...
        'PaperUnits', 'Inches', 'PaperSize', [9, 5], 'renderer','painters')

    ratname = NAcc_ratlist{rat};
    plts = NaN(length(Alignments), 1);
    for a = 1:length(Alignments)
        hh = zeros(length(rewards), 1); %set up to name plots
        for rew = 1:length(rewards)
            da_rat = DAbyrat{rat};
            da_mean = da_rat{rew, a};
            err_rat = DAerrbyrat{rat};
            standarderr = err_rat{rew, a};
            err = [da_mean-standarderr...
                fliplr(da_mean+standarderr)];
            err(isnan(err)) = 0;
            plts(a) = subplot(1, length(Alignments), a);
            hh(rew,1) = plot(T, da_mean, 'Color', mycolors{rew}, ...
                'LineWidth', 0.5); hold on
            fl = fill([T fliplr(T)], err, 'k', 'FaceColor', ...
                mycolors{rew}, 'LineStyle', 'none');
            set(fl, 'facealpha', 0.25);
            set(gca, 'TickDir', 'out'); box off
            xlabel(['Time from ', Alignments{a}]);
            xlim([-1 1.25])
            xticks(-1:0.5:1)
            axis square
        end
        if a == length(Alignments)
            legend(hh, rew_vols, location = 'best')
        end
    end
    sgtitle(ratname)
    allYLim = get(plts, {'YLim'});
    allYLim = cat(2, allYLim{:});
    set(plts, 'YLim', [min(allYLim), max(allYLim)]);
    for a = 1:length(Alignments)
        subplot(1, length(Alignments), a)
        yticks(-2:1:5)
        xline(0, '--k', 'HandleVisibility', 'off')
    end
end

