function plot_fig2e(datadir)
% Plot average event-aligned dopamine release split by delay to reward 
% quartile. Signals are z-scored and baseline corrected before pooling
% across rats (N = 10).
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
files = dir(fullfile(fulldatadir, strcat('*_avgFbyDelay_bc.mat')));
files = {files.name};
fprintf('%i files found\n', length(files));

data_rats = cell(4,5); % delay bin, event
for f=1:length(files)
    load(fullfile(fulldatadir, files{f}), 'data');
    for d=1:4
        for e=1:5
            data_rats{d,e}(f,:) = data{d,e};
        end
    end
end

% drop first event (CPOn)
data_rats = data_rats(:,2:end);

% plotting parameters
mycolor = getColorScheme('delay');
x_range = [-.2 0.8];
T = linspace(-5, 10, 7229);

figure;
set(gcf, 'Color', [1 1 1], 'Units', 'Inches',...
    'Position', [6,4,9.4,2.4],'renderer','painters');
A = {'Offer Cue', 'Reward Port On', ...
    'Reward Cue', 'Reward Delivery'}; % event names
t = tiledlayout(1, length(A));

range = [0 0]; % plot y axis range; for loop below will automatically adjust this
for a=1:length(A)
    nexttile(a);
    for d=1:4
        plotPretty(T, data_rats{d,a}, mycolor{d})
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

