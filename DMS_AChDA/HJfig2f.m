function HJfig2f(datadir)
% Plot event-aligned z-scored DA signal by reward port side relative to
% recording hemisphere (contra/ipsi), averaged across rats (N = 10).
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

folderpath = fullfile(datadir, 'data-published', 'PhotometryData', ...
    'GRAB_DA_DMS');
% load list of DA rats
DArats = load(fullfile(datadir, 'data-published', 'ratlist_DA.mat'), ...
    'ratList');

A = {'CP On', 'Offer Cue', 'Reward Port On', 'Reward Cue', ...
    'Reward', 'Opt Out'};
dataContra = cell(1, length(A));
dataIpsi = dataContra;

% loop through rats to get event-aligned DA signals
for rr=1:length(DArats.ratList)
    ratname = upper(DArats.ratList{rr});
    
    % load event-aligned average z-scored DA split by reward side for each rat
    fname = strcat(ratname, '_avgFbyRewardSide_nbc.mat');
    
    % 'contra' and 'ipsi' are both 6x7229 arrays where rows are task events
    % and columns are timepoints spanning -5 to 10 sec around the events
    load(fullfile(folderpath, fname), 'contra', 'ipsi');
    for a=1:length(A)
        dataContra{a}(rr,:) = contra(a,:);
        dataIpsi{a}(rr,:) = ipsi(a,:);
    end
end

% make a figure
% plotting parameters
xrange = [-0.5, 1]; % show -0.5 to 1 sec around the events
T = linspace(-5, 10, size(dataContra{1},2)); % timepoints
mycolor = {'#ffa052', '#adc178'}; % contra, ipsi
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
    plotPretty(T, dataContra{ind}, mycolor{1});
    plotPretty(T, dataIpsi{ind}, mycolor{2});
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



