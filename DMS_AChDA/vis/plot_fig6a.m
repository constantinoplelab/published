function plot_fig6a(datadir)
% Plot z-scored dopamine and acetylcholine responses on contralateral trials
% aligned to when the side LED turns on, averaged across rats. N = 10 DA
% rats, N = 10 ACh rats.

SON_DA = getSONcontradata(datadir, 'DA');
SON_ACh = getSONcontradata(datadir, 'ACh');

% plotting parameters
mycolor = getColorScheme('side');
mycolor = mycolor.contra;
x_range = [-.1 0.5];
T = linspace(-5, 10, 7229);

figure;
t = tiledlayout(1,2, Padding="compact");
set(gcf, color=[1 1 1], units='Inches',...
    position=[6,4,5,2], renderer='painters');

nexttile(1)
plotPretty(T, SON_DA, mycolor)
ylim([-0.05 0.7])
yticks([0, 0.4])
ylabel('DA')

nexttile(2)
plotPretty(T, SON_ACh, mycolor)
ylim([-0.15 0.47])
yticks([0 0.4])
ylabel('ACh')

for tt=1:2
    nexttile(tt)
    xlim(x_range);
    xline(0, 'k--')
    xlabel('Reward Port On (s)');
    set(gca, 'xtick', [0, 0.2, 0.4]);
    axis square
end

end

function dataC_rats = getSONcontradata(datadir, sensor)

basedir = fullfile(datadir, 'data-published/PhotometryDataAligned');
pdatadir = fullfile(basedir, strcat(sensor, '_DMS'));
fprintf('Loading all rats in %s..\n', pdatadir)
files = dir(fullfile(pdatadir, strcat('*_avgFbySide.mat')));
files = {files.name};
fprintf('%i files found\n', length(files));

dataC_rats = nan(length(files), 7229);
for f=1:length(files)
    load(fullfile(pdatadir, files{f}), 'dataContra');
    dataC_rats(f,:) = dataContra{3}; % side-on is 3rd event
end

end





