function plot_fig6b(datadir)
% Plot latency to peak in dopamine and acetylcholine after the side LED
% turns on (left) and when rats opt out (right) on contralateral trials,
% averaged across rats. N = 10 DA rats, N = 10 ACh rats.
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

region = 'dms';
%% all rats (both separately and simultaneously measured rats)
[ms, ~, ~] = get_avgFlatency(datadir, region, 'da', 'side', 'peak', 0.4);
da_sideon = ms{3}; % DA peak latency at SideOn
da_optout = ms{6}; % DA peak latency at OptOut
[ms, ~, ~] = get_avgFlatency(datadir, region, 'ach', 'side', 'peak', 0.4);
ach_sideon = ms{3}; % ACh peak latency at SideOn
ach_optout = ms{6}; % ACh peak latency at OptOut

figure;
set(gcf, renderer='painters', units='inches', ...
    position = [7,5,3.5,2])
t = tiledlayout(1,2);
nexttile(1);
plot_mean_sem(da_sideon, ach_sideon);
p(1) = ranksum(da_sideon, ach_sideon);
subtitle({sprintf('p=%.3f', p(1)), 'Reward Port On'})

nexttile(2);
plot_mean_sem(da_optout, ach_optout);
p(2) = ranksum(da_optout, ach_optout);
subtitle({sprintf('p=%.3f', p(2)), 'Opt Out'})
ylabel(t, 'Latency (ms)', fontsize=8)

for tt=1:2
    nexttile(tt);
    ylim([160 330])
    yticks([150:50:300])
end

end



function plot_mean_sem(data1, data2)

plot([0.5, 1], [mean(data1), mean(data2)], 'k_', linewidth=1.5);
hold on
errorbar([0.5, 1], [mean(data1), mean(data2)], ...
    [sem(data1), sem(data2)], capsize=0, color='k', linestyle='none', linewidth=1.5)

set(gca, box='off', tickdir='out', xtick=[0.5, 1], xticklabels={'DA', 'ACh'})
    xlim([0 1.5])
end