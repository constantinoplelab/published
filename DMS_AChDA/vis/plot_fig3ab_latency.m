function plot_fig3ab_latency(datadir)
% Plot latency to peak/trough in DA/ACh at offer cue (left) and reward cue (right). 
% N = 10 DA rats, N = 10 ACh rats
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

region = 'dms';
[ms, ~] = run_getAvgFlatency(datadir, region, 'da', 'volume', 'peak');
vol_da = [ms{4,2}; ms{5,2}]; % DA peak latency for 40,80uL at Offer Cue
[ms, ~] = run_getAvgFlatency(datadir, region, 'ach', 'volume', 'trough');
vol_ach = [ms{4,2}; ms{5,2}]; % ACh trough latency for 40,80uL at Offer Cue
[ms, ~] = run_getAvgFlatency(datadir, region, 'da', 'delay', 'peak');
delay_da = cell2mat(ms(:,4)); % DA peak latency at Reward Cue
[ms, ~] = run_getAvgFlatency(datadir, region, 'ach', 'delay', 'trough');
delay_ach = cell2mat(ms(:,4)); % ACh trough latency for longest delay at Reward Cue

figure;
set(gcf, renderer='painters', units='inches', position = [8, 5, 6, 2.5])
t = tiledlayout(1,2);
nexttile(1);
plot_mean_sem(vol_da, vol_ach);
p(1) = ranksum(vol_da, vol_ach);
subtitle({sprintf('DA peak, ACh trough (p=%.4f)', p(1)), ...
    '40,80uL trials at Offer Cue'})

nexttile(2);
plot_mean_sem(delay_da, delay_ach);
p(2) = ranksum(delay_da, delay_ach);
subtitle({sprintf('DA peak, ACh trough (p=%.4f)', p(2)), ...
    'Reward Cue'})

for tt=1:2
    nexttile(tt)
    ylim([150 350])
    yticks(150:100:400)
end
ylabel(t, 'Latency (ms)')

end

function [ms, amp] = ...
    run_getAvgFlatency(datadir, region, sensor, datatype, minmaxtype)

basedir = fullfile(datadir, 'data-published\PhotometryDataAligned');
[region, sensor] = standardizeInputs(region, sensor);
fulldatadir = fullfile(basedir, strcat(sensor, '_', region));

if strcmpi(datatype, 'volume')
    namepattern = '*_avgFbyVol.mat';
elseif strcmpi(datatype, 'delay')
    namepattern = '*_avgFbyDelay.mat';
elseif strcmpi(datatype, 'side')
    namepattern = '*_avgFbySide.mat';
end

files = dir(fullfile(fulldatadir, namepattern));
files = {files.name};
fprintf('%i files found\n', length(files))
ms = {};
amp = {};
tf = 0.5;

for f=1:length(files)
    if ~strcmpi(datatype, 'side')
        load(fullfile(fulldatadir, files{f}), 'data');
    else
        load(fullfile(fulldatadir, files{f}), 'dataContra');
        data = dataContra;
    end
    
    for r=1:size(data,1)
        for c=1:size(data,2)
            [ms{r,c}(f,1), amp{r,c}(f,1)] = ...
                get_avgLatency(data{r,c}, minmaxtype, tf);
        end
    end

end
end

function plot_mean_sem(data1, data2)

plot([mean(data1), mean(data2)], 'k_');
hold on
errorbar([mean(data1), mean(data2)], ...
    [sem(data1), sem(data2)], capsize=0, color='k', linestyle='none')
set(gca, box='off', tickdir='out', xtick=[1 2], xticklabels={'DA', 'ACh'})
    xlim([0.5 2.5])
end
