function plot_fig2f(datadir)
% Plot average event-aligned dopamine release split by reward port side.
% Signals are z-scored and baseline corrected before pooling across rats (N = 10).
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

basedir = fullfile(datadir, 'data-published\PhotometryDataAligned');

region = 'dms';
sensor = 'da';
[region, sensor] = standardizeInputs(region, sensor);
fulldatadir = fullfile(basedir, strcat(sensor, '_', region));

% load data
fprintf('Loading all rats in %s..\n', fulldatadir)
files = dir(fullfile(fulldatadir, strcat('*_avgFbySide_bc.mat')));
files = {files.name};
fprintf('%i files found\n', length(files));

dataC_rats = cell(1,6);
dataI_rats = cell(1,6);
for f=1:length(files)
    load(fullfile(fulldatadir, files{f}), 'dataContra', 'dataIpsi');
    for e=1:6
        dataC_rats{e}(f,:) = dataContra{e};
        dataI_rats{e}(f,:) = dataIpsi{e};
    end
end

% drop first event (CPOn)
dataC_rats = dataC_rats(:,2:end);
dataI_rats = dataI_rats(:,2:end);

% plotting parameters
mycolor = getColorScheme('side');
x_range = [-.2 0.8];
T = linspace(-5, 10, 7229);

figure;
set(gcf, 'Color', [1 1 1], 'Units', 'Inches',...
    'Position', [6,4,10,2.4],'renderer','painters');
A = {'Offer Cue', 'Reward Port On', 'Reward Cue', 'Reward Delivery',...
    'Opt Out'}; % event names
t = tiledlayout(1, length(A));

range = [0 0]; % plot y axis range; for loop below will automatically adjust this
for a=1:length(A)
    nexttile(a);
    plotPretty(T, dataC_rats{a}, mycolor.contra)
    plotPretty(T, dataI_rats{a}, mycolor.ipsi)
    xlim(x_range);
    yl = ylim;
    range = updateYlim(range, yl);
    xlabel(A{a});
    set(gca, 'xtick', [0, 0.3, 0.6]);
    axis square
end
for a = 1:length(A)
    nexttile(a);
    ylim(range);
    xline(0, 'k--')
end
ylabel(t, sprintf('Z-scored %s', sensor))
sgtitle(sprintf('N = %i rats', length(files)))

end

