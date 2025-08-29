function plotFigureS7_TCAsupp(ephysPath, codePath, processedOFCPath)
% Plot TCA error and model, incongruent respones for LED on and LED off

%INPUTS:
%   ephysPath = local path to ephys data downloaded from zenodo 
%   processedEphysPath = local path to saved incongruent response data from
%       running processEphysData_experts

%% Set up paths
s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codePath, 'IgnoreCase', ispc);

if ~onPath % only add code dir to path if it already isn't
    addpath(genpath(codePath))
end

%% general
sem = @(x) std(x,'omitnan')./sqrt(sum(any(~isnan(x), 2)));

setDefaultFigProps

%% TCA path
TCApath = [ephysPath, filesep, 'Expert', filesep, 'OFC', filesep, 'TCA', filesep];

%% plot reconstruction error
error = load([TCApath 'error.mat']);

n_fits = 10;
testComps = 1:20;

figure; hold on
for r = testComps
    x = ones(n_fits,1).*r;
    scatter(x, error.err(:,r), 10, 'k','filled')
end
plot(testComps,mean(error.err),'r')
xticks(0:4:20)
xlabel('# components')
ylabel('Reconstruction error')
set(gca, 'TickDir', 'out'); box off;
set(gcf,'renderer','painter')
set(gcf, 'units', 'centimeters', 'position', [10 10 6 4])


%% plot 8 component TCA model - rank based on error plot

load([TCApath 'model_8comp_timeblockTensor'])

x2 = repmat(-1:.05:0.8, 1, 2);
xvec = [repmat(-1:.05:1, 1, 3) x2];
drawLines = nan(2,5);
drawLines(1,:) = find(xvec==0);
drawLines(2,1:3) = find(xvec == 1);
drawLines(2,4:5) = find(xvec == 0.8, 2, 'last');

titles = {'Neuron' 'Time' 'Block'};

plotTCA_ephysFig(model, drawLines,  ...
    'Plottype', {'bar', 'line', 'line'}, ...
    'Modetitles', titles)
set(gcf, 'units', 'centimeters', 'position', [10 10 9 8])
set(gcf,'renderer','painter')

%% plot incongruent responses for each cluster for LED on and LED off alignments

timeBin = [-4 8]; % time window around each event saved in the SU structs
x = timeBin(1):0.05:timeBin(2);
wndw = [-0.5 1];
w = arrayfun(@(y) find(x == wndw(y)), 1:2); %indices of window for computing significance
xvec = w(1):w(2);

%load incongruent responses for expert OFC
incongruentResponses = load([processedOFCPath, 'incongruentResponses.mat']);

% load cluster info for each cell
tcaInfo_expert = load([TCApath 'info_8comp.mat']);
tcaInfo_expert = tcaInfo_expert.info;

clusts = unique(tcaInfo_expert.TCAclust);
nc = length(clusts); 

%make combo table to pull cell responses by cluster
expertTable = join(incongruentResponses.info, tcaInfo_expert);


figure; hold on
tiledlayout(nc, 4, "TileSpacing", 'compact')

labels = {'Trial start' 'LED on' 'LED off' 'Reward', 'Opt out'};
for c = 1:nc  

    % now plot cluster incongruent responses for trial start and reward
    for  a = [2 3]

        clustPre_pref = incongruentResponses.pre_pref{a}(expertTable.TCAclust == c-1, :);
        clustPre_pref_Avg = mean(clustPre_pref, 1, 'omitnan');
        clustPre_pref_Sem = sem(clustPre_pref);

        clustPre_nonpref = incongruentResponses.pre_nonpref{a}(expertTable.TCAclust == c-1, :);
        clustPre_nonpref_Avg = mean(clustPre_nonpref, 1, 'omitnan');
        clustPre_nonpref_Sem = sem(clustPre_nonpref);

        clustPost_pref = incongruentResponses.post_pref{a}(expertTable.TCAclust == c-1, :);
        clustPost_pref_Avg = mean(clustPost_pref, 1, 'omitnan');
        clustPost_pref_Sem = sem(clustPost_pref);

        clustPost_nonpref = incongruentResponses.post_nonpref{a}(expertTable.TCAclust == c-1, :);
        clustPost_nonpref_Avg = mean(clustPost_nonpref, 1, 'omitnan');
        clustPost_nonpref_Sem = sem(clustPost_nonpref);

        % run permutation test on each time bin
        pExp_post = arrayfun(@(y) permute_test(clustPost_nonpref(:,y), ...
            clustPost_pref(:,y), 1000, 1), xvec);
        pExp_pre = arrayfun(@(y) permute_test(clustPre_nonpref(:,y), ...
            clustPre_pref(:,y), 1000, 1), xvec);

        % find significant bins
        pExp_pre_sig = pExp_pre < 0.05;
        pExp_post_sig = pExp_post < 0.05;


        %find 3 consecutive significant bins
        sig_pre = false(1, length(xvec));
        a0 = [0 pExp_pre_sig 0]; % add 0 at beginning and end to find pattern. also shifts index so start is the first true value
        startTrue = strfind(a0,[0 1]);
        endTrue = strfind(a0,[1 0]);
        numSigBins = endTrue - startTrue;
        use = find(numSigBins > 2);
        good = arrayfun(@(x) startTrue(use(x)):endTrue(use(x))-1, 1:length(use), ...
            'uniformoutput', false);
        sig_pre(cell2mat(good)) = true;


        %plot
        ax = nexttile;
        shadedErrorBar(x, clustPre_pref_Avg, clustPre_pref_Sem, 'lineprops', ...
            {'color', [0.65 0.65 0.65], 'linewidth', 1})
        hold on
        shadedErrorBar(x, clustPre_nonpref_Avg, clustPre_nonpref_Sem, 'lineprops', ...
            {'color', 'k', 'linewidth', 1})
        xlim([-0.5 1])
        % ylim([-0.5 1])
        yl = ylim;
        scatter(x(xvec(sig_pre)), yl(2), '*k')
        plot([0 0], [yl(1) yl(2)], '--k', 'linewidth', 0.25)
        hold on
        % set(gca, 'ytick', [])
        if c == nc
            set(gca, 'xtick', [-0.5, 0, 1])
            xlabel(labels{a})
        else
            set(gca, 'xtick', [])
        end
        set(gca, 'TickDir', 'out'); box off;
        ax.YRuler.TickLabelGapOffset = 1;
    

        %find 3 consecutive significant bins
        sig_post = false(1, length(xvec));
        a0 = [0 pExp_post_sig 0]; % add 0 at beginning and end to find pattern. also shifts index so start is the first true value
        startTrue = strfind(a0,[0 1]);
        endTrue = strfind(a0,[1 0]);
        numSigBins = endTrue - startTrue;
        use = find(numSigBins > 2);
        good = arrayfun(@(x) startTrue(use(x)):endTrue(use(x))-1, 1:length(use), ...
            'uniformoutput', false);
        sig_post(cell2mat(good)) = true;


        ax2 = nexttile;
        shadedErrorBar(x, clustPost_pref_Avg, clustPost_pref_Sem, 'lineprops', ...
            {'color', [0.65 0.65 0.65], 'linewidth', 1})
        hold on
        shadedErrorBar(x, clustPost_nonpref_Avg, clustPost_nonpref_Sem, 'lineprops', ...
            {'color', 'k', 'linewidth', 1})
        xlim([-0.5 1])
        ylim([yl])
        scatter(x(xvec(sig_post)), yl(2), '*k')
        plot([0 0], [yl(1) yl(2)], '--k', 'linewidth', 0.25)
        hold on
        set(gca, 'ytick', [])
        if c == nc
            set(gca, 'xtick', [-0.5, 0, 1])
            xlabel(labels{a})
        else
            set(gca, 'xtick', [])
        end
        set(gca, 'TickDir', 'out'); box off;
        ax2.YRuler.TickLabelGapOffset = 1;

    end
    
end
set(gcf,'renderer','painter')
set(gcf, 'units', 'centimeters', 'position', [10 10 10 12])


