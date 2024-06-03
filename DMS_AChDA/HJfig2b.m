function HJfig2b(datadir)
% Plot event-aligned z-scored DA signal by volume in mixed blocks,
% averaged across rats (N = 10).
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

folderpath = fullfile(datadir, 'data-published', 'PhotometryData', ...
    'GRAB_DA_DMS');
% load list of DA rats
DArats = load(fullfile(datadir, 'data-published', 'ratlist_DA.mat'), ...
    'ratList');

rewards = [5,10,20,40,80];
A = {'CP On', 'Offer Cue', 'Reward Port On', 'Reward Cue', ...
    'Reward', 'Opt Out'};
dataAllRats = cell(length(rewards), length(A));

% loop through rats to get event-aligned DA signals
for rr=1:length(DArats.ratList)
    ratname = upper(DArats.ratList{rr});
    
    % load event-aligned average z-scored DA split by volume for each rat
    fname = strcat(ratname, '_avgFbyVolMixedBlock_nbc.mat');

    % 1x5 cell array where each cell corresponds to different reward volumes
    % in ascending order.
    % Each cell array is 6x7229 where rows are different task events (CP
    % On, Offer Cue, Reward Port On, Reward Cue, Reward, Opt Out) and
    % columns are timepoints spanning -5 to 10 sec around the events.
    load(fullfile(folderpath, fname), 'data');
    for rew=1:length(rewards)
        for a=1:length(A)
            dataAllRats{rew,a}(rr,:) = data{rew}(a,:);
        end
    end
end

% make a figure
% plotting parameters
xrange = [-0.5, 1]; % show -0.5 to 1 sec around the events
T = linspace(-5, 10, size(dataAllRats{1,1},2)); % timepoints
mycolor = {'#bb7af7', '#9a5bd0', '#793cab', '#591e86', '#390061'}; 
% task events to visualize
Avis = {'Offer Cue', 'Reward Port On', 'Reward Cue', 'Opt Out'};

figure;
t = tiledlayout(1,4, Padding="compact");
set(gcf, units='inches', position=[6,5,8.8,2.65])
yrange = [0 0]; % code below will adjust this to make y-axis same across all plots
for aa=1:length(Avis)
    nexttile(aa);
    hold on
    ind = find(cellfun(@(x) strcmpi(x, Avis{aa}), A));
    for rew=1:length(rewards)
        plotPretty(T, dataAllRats{rew,ind}, mycolor{rew});
    end
    xlim(xrange);
    xlabel(sprintf('%s (s)', Avis{aa}));
    xline(0, 'k--')
    yl = ylim;
    yrange = updateYlim(yrange, yl);
    axis square
end

for aa=1:length(Avis)
    nexttile(aa);
    ylim(yrange);
end

ylabel(t, 'Z-scored DA', fontsize=13)


end