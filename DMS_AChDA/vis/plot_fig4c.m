function plot_fig4c(datadir)
% Plot average time to opt-out in mixed blocks for rats tested with catch
% probabilities of 0.15 and 0.25. N = 53 rats.
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

% Load wait time data of 53 rats 
fpath = fullfile(datadir, 'data-published/WTbyPcatch.mat');
load(fpath, 'WT')

figure;
set(gcf, units='inches', renderer='painters',...
    position=[7,4,1.6,2])
colors = {'#BCBEC0', '#231F20'};

wt_avg(:,1) = cellfun(@(x) mean(x(:,1), 'omitnan'), WT);
wt_avg(:,2) = cellfun(@(x) mean(x(:,2), 'omitnan'), WT);
wt_sem(:,1) = cellfun(@(x) sem(x(:,1), 'omitnan'), WT);
wt_sem(:,2) = cellfun(@(x) sem(x(:,2), 'omitnan'), WT);

for pp=1:2
    plot(1:5, wt_avg(:,pp), color=colors{pp}, linewidth=1)
    hold on
    errorbar(1:5, wt_avg(:,pp), wt_sem(:,pp), color=colors{pp}, CapSize=0, ...
        LineWidth=1)
end

set(gca, box='off', tickdir='out', fontsize=8, xtick=1:5, ...
    xticklabels={'5', '10', '20', '40', '80'})
xlim([0.5 5.5])
ylim([9 17])
xlabel('Offered reward')
ylabel('Time to opt-out (s)')
subtitle(sprintf('N = %i rats', length(WT{1})))

% stats
p_signrank = nan(1,5);
for rew=1:5
    p_signrank(rew) = signrank(WT{rew}(:,1), WT{rew}(:,2));
end
disp(p_signrank)



end