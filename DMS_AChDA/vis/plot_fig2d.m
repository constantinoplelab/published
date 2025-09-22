function plot_fig2d(datadir)
% Plot average event-aligned dopamine release during mixed blocks for
% different reward offers. Signals are z-scored and baseline corrected
% before pooling across rats (N = 10).
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
files = dir(fullfile(fulldatadir, strcat('*_avgFbyVol_bc.mat')));

files = {files.name};
fprintf('%i files found\n', length(files));

data_rats = cell(5,6); % reward, event
for f=1:length(files)
    load(fullfile(fulldatadir, files{f}), 'data');
    for rew=1:5
        for e=1:6
            data_rats{rew,e}(f,:) = data{rew,e};
        end
    end
end

% drop first event (CPOn)
data_rats = data_rats(:,2:end);

% plotting parameters
mycolor = getColorScheme('volume');
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
    for rew=1:5
        plotPretty(T, data_rats{rew,a}, mycolor{rew})
    end
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

