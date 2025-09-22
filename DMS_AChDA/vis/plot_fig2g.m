function plot_fig2g(datadir)
% Plot event-aligned ACh responses (baseline corrected), split by relevant
% conditions at each event, averaged across rats (N = 10).
% 1. Offer Cue: volume
% 2. Reward Port On: side
% 3. Reward Cue: delay
% 4. Reward Delivery: delay
% 5. Opt Out: side
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

ratlist = loadRats(datadir, 'ach');

[ACh_offercue, ACh_sideon, ACh_rewardcue, ACh_reward, ACh_optout] = ...
    get_avgFl_allevents(datadir, ratlist, 'ACh', 'DMS');
plotavgfl(ACh_offercue, ACh_sideon, ...
    ACh_rewardcue, ACh_reward, ACh_optout)

nexttile(1); ylim([-0.25 0.15]); yticks([-0.2 0])
nexttile(2); yticks([0, 0.4])
nexttile(3); ylim([-0.3 0.15]); yticks([-0.2, 0])
nexttile(4); ylim([-0.3 0.15]); yticks([-0.2, 0])
nexttile(5); ylim([-0.05 0.15]); yticks([0, 0.1])
end

function plotavgfl(avgf_offercue, avgf_sideon, ...
    avgf_sideoff, avgf_reward, avgf_optout)

%%
T = linspace(-5, 10, 7229);
trange = [-0.2, 0.8];
figure
t = tiledlayout(1,5, Padding="compact");
set(gcf, units='inches', renderer='painters', ...
    position=[7,3,7.6,1.5])

%% 1. Offer Cue by volume
mycolors = getColorScheme('volume');
rewards = [5, 10, 20, 40, 80];

nexttile(1)
hold on
for rew=1:length(rewards)
    plotPretty(T, avgf_offercue{rew}, mycolors{rew});
end

%% 2. Reward Port On and Opt Out by side
mycolors = getColorScheme('side');

nexttile(2)
hold on
plotPretty(T, avgf_sideon.contra, mycolors.contra);
plotPretty(T, avgf_sideon.ipsi, mycolors.ipsi);

nexttile(5)
hold on
plotPretty(T, avgf_optout.contra, mycolors.contra);
plotPretty(T, avgf_optout.ipsi, mycolors.ipsi);

%% 3. Reward Cue and Reward Delivery by delay
mycolors = getColorScheme('delay');
nbin=4;

for bb=1:nbin
    nexttile(3)
    hold on
    plotPretty(T, avgf_sideoff{bb}, mycolors{bb});
    
    nexttile(4)
    hold on
    plotPretty(T, avgf_reward{bb}, mycolors{bb});
end

%%
for tt=1:5
    nexttile(tt)
    xlim(trange)
    xline(0, 'k--')
    set(gca, fontsize=8, xtick=[0,0.3,0.6])
end

xlabel(t, 'Time from event (s)')
ylabel(t, 'Z-scored F')
drawnow

end



