function plotFigure3_expertEphys(ephysPath, processedEphysPath, codePath)
% Plot figure 3. YOU MUST RUN processEphysData_experts ON THE OFC
% ELECTROPHYSIOLOGY DATA BEFORE RUNNING THIS FUNCTION. Plots average firing
% rates for 20ul trials before and after the first incongruent trial.

% INPUTS:
%   ephysPath = local path to .mat files for ephys data downloaded from
%       zenodo
%   processedEphysPath = local path to saved output from
%       processEphysData_experts using OFC data files
%   codePath = local path to code saved from github


%% Set up paths
s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codePath, 'IgnoreCase', ispc);

if ~onPath % only add code dir to path if it already isn't
    addpath(genpath(codePath))
end

%% general

timeBin = [-4 8]; % time window around each event saved in the SU structs
x = timeBin(1):0.05:timeBin(2);
wndw = [-0.5 1];
w = arrayfun(@(y) find(x == wndw(y)), 1:2); %indices of window for computing significance
xvec = w(1):w(2);

sem = @(x) std(x,'omitnan')./sqrt(sum(any(~isnan(x), 2)));

setDefaultFigProps

%% OFC path
OFCpath = [ephysPath, filesep, 'Expert', filesep, 'OFC', filesep];

%% load incongruent responses
incongruentResponses = load([processedEphysPath, 'incongruentResponses.mat']);

%% averages
nAlign = 4;

%averages
prePref_Avg = arrayfun(@(x) mean(incongruentResponses.pre_pref{x}, 1, 'omitnan'), ...
    1:nAlign, 'UniformOutput', false);
prePref_Sem = arrayfun(@(x) sem(incongruentResponses.pre_pref{x}), 1:nAlign, ...
    'uniformoutput', false);

preNonpref_Avg = arrayfun(@(x) mean(incongruentResponses.pre_nonpref{x}, 1, 'omitnan'), ...
    1:nAlign, 'UniformOutput', false);
preNonpref_Sem = arrayfun(@(x) sem(incongruentResponses.pre_nonpref{x}), 1:nAlign, ...
    'uniformoutput', false);

postPref_Avg = arrayfun(@(x) mean(incongruentResponses.post_pref{x}, 1, 'omitnan'), ...
    1:nAlign, 'UniformOutput', false);
postPref_Sem = arrayfun(@(x) sem(incongruentResponses.post_pref{x}), 1:nAlign, ...
    'uniformoutput', false);

postNonpref_Avg = arrayfun(@(x) mean(incongruentResponses.post_nonpref{x}, 1, 'omitnan'), ...
    1:nAlign, 'UniformOutput', false);
postNonpref_Sem = arrayfun(@(x) sem(incongruentResponses.post_nonpref{x}), 1:nAlign, ...
    'uniformoutput', false);


%% compute p-values for each time bin

% pExp_post = arrayfun(@(x) permute_test(mean(incongruentResponses.post_nonpref{x}...
%     (:, w(1):w(2)), 2, 'omitnan'), mean(incongruentResponses.post_pref{x}...
%     (:, w(1):w(2)), 2, 'omitnan'), 10000), 1:nAlign);
% pExp_pre = arrayfun(@(x) permute_test(mean(incongruentResponses.pre_nonpref{x}...
%     (:, w(1):w(2)), 2, 'omitnan'), mean(incongruentResponses.pre_pref{x}...
%     (:, w(1):w(2)), 2, 'omitnan'), 10000), 1:nAlign);

pExp_post = cell(1, nAlign);
pExp_pre = pExp_post;
for jj = 1:nAlign
    pExp_post{jj} = arrayfun(@(y) permute_test(...
        incongruentResponses.post_nonpref{jj}(:,y), ...
        incongruentResponses.post_pref{jj}(:,y), 1000, 1), xvec);
    pExp_pre{jj} = arrayfun(@(y) permute_test(...
        incongruentResponses.pre_nonpref{jj}(:,y), ...
        incongruentResponses.pre_pref{jj}(:,y), 1000, 1), xvec);
end

% find significant bins
pExp_pre_sig = arrayfun(@(x) pExp_pre{x} < 0.05, 1:nAlign, ...
    'uniformoutput', false);
pExp_post_sig = arrayfun(@(x) pExp_post{x} < 0.05, 1:nAlign, ...
    'uniformoutput', false);


%% example block significant cell
filename = 'S021_2021-09-22.mat';
cellID = 150;

load(strcat(OFCpath, filename), 'S', 'SU');
ids = extractfield(SU, 'cluster_id');
use = ids == cellID;

hmat = SU(use).hmat.Rew;
byBlock = arrayfun(@(x) mean(hmat(S.Block == x, :), 1, 'omitnan'), ...
    [2, 3], 'UniformOutput', false);
byBlock_sem = arrayfun(@(x) sem(hmat(S.Block == x, :)), [2 3], ...
    'UniformOutput', false);

%% plot incongruent responses

figure; hold on
tiledlayout(2, nAlign+1, 'TileSpacing', 'compact')

%plot example cell with significantly different firing rates between low
%and high blocks
nexttile
shadedErrorBar(x, byBlock{1}, byBlock_sem{1}, 'lineprops', {'color', ...
    'r', 'linewidth', 1})
shadedErrorBar(x, byBlock{2}, byBlock_sem{2}, 'lineprops', {'color', ...
    'b', 'linewidth', 1})
set(gca, 'TickDir', 'out'); box off;
xlim([-0.5 2]);
xlabel('Time from reward (s)')
ylabel('Firing rate (Hz)')
axis square
ax1 = gca;
ax1.YRuler.TickLabelGapOffset = 1;
title('\rm Example cell')

%pre-incongruent
labels = {'Trial start' 'LED on' 'LED off' 'Reward'};
for ii = 1:nAlign

    sig = false(1, length(xvec));
    %find 3 consecutive significant bins
    a0 = [0 pExp_pre_sig{ii} 0]; % add 0 at beginning and end to find pattern. also shifts index so start is the first true value
    startTrue = strfind(a0,[0 1]);
    endTrue = strfind(a0,[1 0]);
    numSigBins = endTrue - startTrue;
    use = find(numSigBins > 2);
    good = arrayfun(@(x) startTrue(use(x)):endTrue(use(x))-1, 1:length(use), ...
        'uniformoutput', false);
    % sig_pre = xvec(cell2mat(good));
    sig(cell2mat(good)) = true;
  
    % plot
    ax(ii) = nexttile;
    shadedErrorBar(x(xvec), prePref_Avg{ii}(xvec), prePref_Sem{ii}(xvec), ...
        'lineprops', {'color', [0.65 0.65 0.65], 'linewidth', 1})
    hold on
    shadedErrorBar(x(xvec), preNonpref_Avg{ii}(xvec), preNonpref_Sem{ii}(xvec), ...
        'lineprops', {'color', 'k', 'linewidth', 1})
    scatter(x(xvec(sig)), 0.4, '*k')
    xline(0, '--', 'linewidth', 0.5, 'color', [0.7 0.7 0.7])
    set(gca, 'TickDir', 'out'); box off;
    axis square
    ax(ii).YRuler.TickLabelGapOffset = 1;
    title(['\rm' labels{ii}])
end
ylabel(ax(1), 'Firing rate (Hz)')
linkaxes(ax)
xlim([-0.5 1])
ylim([-0.2 0.5])

%blank tile for schematic
nexttile
axis square


for ii = 1:nAlign
    
    sig = false(1, length(xvec));
    %find 3 consecutive significant bins
    a0 = [0 pExp_post_sig{ii} 0]; % add 0 at beginning and end to find pattern. also shifts index so start is the first true value
    startTrue = strfind(a0,[0 1]);
    endTrue = strfind(a0,[1 0]);
    numSigBins = endTrue - startTrue;
    use = find(numSigBins > 2);
    good = arrayfun(@(x) startTrue(use(x)):endTrue(use(x))-1, 1:length(use), ...
        'uniformoutput', false);
    sig(cell2mat(good)) = true;
    % sig_post = xvec(cell2mat(good));

    ax(ii) = nexttile;
    shadedErrorBar(x(xvec), postPref_Avg{ii}(xvec), postPref_Sem{ii}(xvec), ...
        'lineprops', {'color', [0.65 0.65 0.65], 'linewidth', 1})
    hold on
    shadedErrorBar(x(xvec), postNonpref_Avg{ii}(xvec), postNonpref_Sem{ii}(xvec), ...
        'lineprops', {'color', 'k', 'linewidth', 1})
    scatter(x(xvec(sig)), 0.4, '*k')
    xline(0, '--', 'linewidth', 0.5, 'color', [0.7 0.7 0.7])
    set(gca, 'TickDir', 'out'); box off;
    axis square
    ax(ii).YRuler.TickLabelGapOffset = 1;
    title(['\rm' labels{ii}])
end
ylabel(ax(1), 'Firing rate (Hz)')
linkaxes(ax)
xlim([-0.5 1])
ylim([-0.2 0.5])

set(gcf, 'units', 'centimeters', 'position', [10 10 20 8])
set(gcf,'renderer','painter')
