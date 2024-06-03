function HJfig2h(datadir)
% Plot z-scored ACh signal split by task variable most strongly encoded
% by DA at each task event (e.g., by volume at Offer Cue and reward side at
% Reward Port On), averaged across rats (N = 7).
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

folderpath = fullfile(datadir, 'data-published', 'PhotometryData', ...
    'GRAB_ACh_DMS');
% load list of ACh rats
AChrats = load(fullfile(datadir, 'data-published', 'ratlist_ACh.mat'), ...
    'ratList');


% task events to visualize
Avis = {'Offer Cue', 'Reward Port On', 'Reward Cue', 'Opt Out'};

dataAllRats = cell(1, length(Avis));

% loop through rats to get event-aligned ACh signals
for rr=1:length(AChrats.ratList)
    ratname = upper(AChrats.ratList{rr});
    
    % Load data files relevant for each event

    % Offer Cue: load data split by volume
    rewards = [5, 10, 20, 40, 80];
    fname = strcat(ratname, '_avgFbyVolMixedBlock_nbc.mat');
    load(fullfile(folderpath, fname), 'data');
    ind = cellfun(@(x) strcmpi(x, 'Offer Cue'), Avis);
    for rew=1:length(rewards)
        dataAllRats{ind}{rew}(rr,:) = data{rew}(2,:); 
    end

    % Reward Port On, Opt Out: load data split by reward port side
    fname = strcat(ratname, '_avgFbyRewardSide_nbc.mat');
    load(fullfile(folderpath, fname), 'contra', 'ipsi');
    ind = cellfun(@(x) strcmpi(x, 'Reward Port On'), Avis);
    dataAllRats{ind}{1}(rr,:) = contra(3,:);
    dataAllRats{ind}{2}(rr,:) = ipsi(3,:);
    ind = cellfun(@(x) strcmpi(x, 'Opt Out'), Avis);
    dataAllRats{ind}{1}(rr,:) = contra(6,:);
    dataAllRats{ind}{2}(rr,:) = ipsi(6,:);

    % Reward Cue: load data split by reward delay quartile
    nbins = 4;
    fname = strcat(ratname, '_avgFbyRewardDelay_nbc.mat');
    load(fullfile(folderpath, fname), 'data');
    ind = cellfun(@(x) strcmpi(x, 'Reward Cue'), Avis);
    for nn=1:nbins
        dataAllRats{ind}{nn}(rr,:) = data{nn}(4,:);
    end
end

% make a figure
% plotting parameters
xrange = [-0.5, 1]; % show -0.5 to 1 sec around the events
T = linspace(-5, 10, size(dataAllRats{1}{1},2)); % timepoints

figure;
t = tiledlayout(1,4, Padding="compact");
set(gcf, units='inches', position=[6,5,8.8,2.65])
mycolors{1} = {'#bb7af7', '#9a5bd0', '#793cab', '#591e86', '#390061'};
mycolors{2} = {'#ffa052', '#adc178'};
mycolors{3} = {'#adadad', '#8d8d8d', '#6e6e6e', '#4f4f4f'};
mycolors{4} = mycolors{2};

for aa=1:length(Avis)
    nexttile(aa);
    hold on
    for jj=1:length(dataAllRats{aa})
        plotPretty(T, dataAllRats{aa}{jj}, mycolors{aa}{jj});
    end
    xlim(xrange);
    xlabel(sprintf('%s (s)', Avis{aa}));
    xline(0, 'k--')
    if aa==3
        ylim([-0.4 0.4])
    else
        yl = ylim;
        ylim([min(yl)-0.005 max(yl)+0.005]);
    end
    axis square
end

end




