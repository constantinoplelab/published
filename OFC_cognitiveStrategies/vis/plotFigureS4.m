function plotFigureS4(ephysPath)

%paths
expertPath = [ephysPath, filesep, 'Expert', filesep];
naivePath = [ephysPath, filesep, 'Naive', filesep];

error_expert = load([expertPath, 'TCA', filesep, 'error']);
error_naive = load([naivePath, 'TCA', filesep, 'error']);

n_fits = 10;
test = 1:20;

figure; hold on
tiledlayout(1, 2, 'TileSpacing', 'compact')

nexttile
for r = 1:length(test)
    x = ones(n_fits,1).*test(r);
    scatter(x,error_expert.err(:,r), 10, 'k','filled')
    hold on
end
plot(test,mean(error_expert.err),'r')
xticks(0:4:20)
xlabel('# components')
ylabel('Reconstruction error')
title('\rm Experts')

nexttile
for r = 1:length(test)
    x = ones(n_fits,1).*test(r);
    scatter(x,error_naive.err(:,r), 10, 'k','filled')
    hold on
end
plot(test,mean(error_naive.err),'r')
xticks(0:4:20)
xlabel('# components')
ylabel('Reconstruction error')
title('\rm Naive')
set(gcf, 'renderer', 'painters')
set(gcf, 'units', 'centimeters', 'position', [10 10 16 4])

% plot TCA model - experts
load([expertPath, 'TCA', filesep, 'model_8comp_timeblockTensor']);
load([expertPath, 'TCA', filesep, 'info_8comp']);

x2 = repmat(-1:.05:0.8, 1, 2);
xvec = [repmat(-1:.05:1, 1, 3) x2];
drawLines = nan(2,5);
drawLines(1,:) = find(xvec==0);
drawLines(2,1:3) = find(xvec == 1);
drawLines(2,4:5) = find(xvec == 0.8, 2, 'last');

titles = {'Neuron' 'Time' 'Block'};

plotTCAmodels(model, drawLines,  ...
    'Plottype', {'bar', 'line', 'line'}, ...
    'Modetitles', titles)
sgtitle('\rm Experts')
set(gcf, 'renderer', 'painters')
set(gcf, 'units', 'centimeters', 'position', [10 10 8.5 8])

%plot TCA model - naive
load([naivePath, 'TCA', filesep, 'model_8comp_timeblockTensor']);
load([naivePath, 'TCA', filesep, 'info_8comp']);

x2 = repmat(-1:.05:0.8, 1, 2);
xvec = [repmat(-1:.05:1, 1, 3) x2];
drawLines = nan(2,5);
drawLines(1,:) = find(xvec==0);
drawLines(2,1:3) = find(xvec == 1);
drawLines(2,4:5) = find(xvec == 0.8, 2, 'last');

titles = {'Neuron' 'Time' 'Block'};

plotTCAmodels(model, drawLines,  ...
    'Plottype', {'bar', 'line', 'line'}, ...
    'Modetitles', titles)
sgtitle('\rm Naive')
set(gcf, 'renderer', 'painters')
set(gcf, 'units', 'centimeters', 'position', [10 10 8.5 8])

