function plotFigure7_naiveEphys(ephysPath, codePath)
% Plot figure 7. Naive behavior for all recording sessions, ATI comparison
% between experts and naive, responses to 20ul trials before and after the
% first incongruent trial for low and high ATI sessions

% INPUTS:
%   ephysPath = local path to .mat files for ephys data downloaded from
%       zenodo
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

%% wait time curves for naive ephys sessions
behPath = [ephysPath, filesep, 'Naive',  filesep, 'Behavior', filesep];
load([behPath 'ratList_naive.mat']); %list of naive rats (ephys and behavior only)
nrats = length(ratList);

for rr = 1:nrats
    
    N = load([behPath 'ratTrial_', ratList{rr}, '.mat']);

    % wait time curves -- rat data
    [hi(rr), lo(rr), mix(rr)] = wtcurves_SS(N.A);

end

naive.hi = cell2mat(arrayfun(@(x) hi(x).wt, 1:nrats, 'uniformoutput', false)');
naive.lo = cell2mat(arrayfun(@(x) lo(x).wt, 1:nrats, 'uniformoutput', false)');
naive.mix =  cell2mat(arrayfun(@(x) mix(x).wt, 1:nrats, 'uniformoutput', false)');

avgNm = mean(naive.mix, 'omitnan');
semNm = sem(naive.mix);

avgNh = mean(naive.hi, 'omitnan');
semNh = sem(naive.hi);

avgNl = mean(naive.lo, 'omitnan');
semNl = sem(naive.lo);

%kruskal wallis test for 20ul block effect on 20ul wait times
p_wt20 = kruskalwallis([naive.mix(:,3), naive.hi(:,3), naive.lo(:,3)], [], 'off');

%wilcoxon signrank for other volumes
p_wt5 = signrank(naive.mix(:,1), naive.lo(:,1));
p_wt10 = signrank(naive.mix(:,2), naive.lo(:,2));
p_wt40 = signrank(naive.mix(:,4), naive.hi(:,4));
p_wt80 = signrank(naive.mix(:,5), naive.hi(:,5));

%% ATI
prc = 40; %percentile to split blockiness

naivePath = [ephysPath, filesep, 'Naive', filesep, 'OFC'];
tablePath = [ephysPath, filesep, 'ATI.csv'];

%incongruent responses for block naive rats split by ATI
[low, high] = incongruentResponse_twentyPreVsPost_byATI(naivePath, ...
    tablePath, prc);

T_combined = readtable(tablePath);

ATI = T_combined.ATI;
group = string(T_combined.Group);

naive_ATI = ATI(strcmp(group, 'Naive'));
expert_ATI = ATI(strcmp(group, 'Expert'));
p_ATI = ranksum(expert_ATI, naive_ATI);


%% plot 

x_exp = linspace(min(expert_ATI), max(expert_ATI));
x_naive = linspace(min(naive_ATI), max(naive_ATI));

figure; hold on
tiledlayout(1, 4, 'TileSpacing', 'compact')

%wait time curves
nexttile
shadedErrorBar(1:5, avgNm, semNm, 'lineprops', {'k', 'linewidth', 1})
hold on
shadedErrorBar(1:5, avgNh, semNh, 'lineprops', {'-r', 'linewidth', 1})
shadedErrorBar(1:5, avgNl, semNl, 'lineprops', {'-b', 'linewidth', 1})
set(gca, 'xTick', 1:5);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
ylim([10 16])
xlabel('Reward offer')
ylabel('Mean wait time (s)')
ax1 = gca;
ax1.YRuler.TickLabelGapOffset = 1;
axis square

%z-latent aligned to incongruent -- DPR code
nexttile

%ATI demo -- DPR code
nexttile

ax = nexttile;
histogram(expert_ATI, 'facecolor', 'k')
xline(median(expert_ATI), '--k', 'LineWidth', 1)
hold on
histogram(naive_ATI, 'facecolor', 'm')
xline(median(naive_ATI), '--m', 'LineWidth', 1)
set(gca, 'TickDir', 'out'); box off;
axis square
ax.YRuler.TickLabelGapOffset = 1;
ylabel('N sessions')
xlabel('Across-trial influence')

set(gcf, 'units', 'centimeters', 'position', [10 10 20 4])
set(gcf,'renderer','painter')


%% plot incongruent responses - low ATI

nAlign = 4;

%averages
prePref_Avg = arrayfun(@(x) mean(low.pre_pref{x}, 1, 'omitnan'), ...
    1:nAlign, 'UniformOutput', false);
prePref_Sem = arrayfun(@(x) sem(low.pre_pref{x}), 1:nAlign, ...
    'uniformoutput', false);

preNonpref_Avg = arrayfun(@(x) mean(low.pre_nonpref{x}, 1, 'omitnan'), ...
    1:nAlign, 'UniformOutput', false);
preNonpref_Sem = arrayfun(@(x) sem(low.pre_nonpref{x}), 1:nAlign, ...
    'uniformoutput', false);

postPref_Avg = arrayfun(@(x) mean(low.post_pref{x}, 1, 'omitnan'), ...
    1:nAlign, 'UniformOutput', false);
postPref_Sem = arrayfun(@(x) sem(low.post_pref{x}), 1:nAlign, ...
    'uniformoutput', false);

postNonpref_Avg = arrayfun(@(x) mean(low.post_nonpref{x}, 1, 'omitnan'), ...
    1:nAlign, 'UniformOutput', false);
postNonpref_Sem = arrayfun(@(x) sem(low.post_nonpref{x}), 1:nAlign, ...
    'uniformoutput', false);

%stats
pExp_post = cell(1, nAlign);
pExp_pre = pExp_post;
for jj = 1:nAlign
    pExp_post{jj} = arrayfun(@(y) permute_test(low.post_nonpref{jj}(:,y), ...
        low.post_pref{jj}(:,y), 1000, 1), xvec);
    pExp_pre{jj} = arrayfun(@(y) permute_test(low.pre_nonpref{jj}(:,y), ...
        low.pre_pref{jj}(:,y), 1000, 1), xvec);
end

% find significant bins
pExp_pre_sig = arrayfun(@(x) pExp_pre{x} < 0.05, 1:nAlign, ...
    'uniformoutput', false);
pExp_post_sig = arrayfun(@(x) pExp_post{x} < 0.05, 1:nAlign, ...
    'uniformoutput', false);

%plot
figure; hold on
tiledlayout(2, nAlign, 'TileSpacing', 'compact')

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
    scatter(x(xvec(sig)), 1, '*k')
    xline(0, '--', 'linewidth', 0.5, 'color', [0.7 0.7 0.7])
    set(gca, 'TickDir', 'out'); box off;
    axis square
    ax(ii).YRuler.TickLabelGapOffset = 1;
    title(['\rm' labels{ii}])
end
ylabel(ax(1), 'Firing rate (Hz)')
linkaxes(ax)
xlim([-0.5 1])
ylim([-0.5 1])

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
    scatter(x(xvec(sig)), 1, '*k')
    xline(0, '--', 'linewidth', 0.5, 'color', [0.7 0.7 0.7])
    set(gca, 'TickDir', 'out'); box off;
    axis square
    ax(ii).YRuler.TickLabelGapOffset = 1;
    title(['\rm' labels{ii}])
end
ylabel(ax(1), 'Firing rate (Hz)')
linkaxes(ax)
xlim([-0.5 1])
ylim([-0.4 1])

set(gcf, 'units', 'centimeters', 'position', [10 10 16 8])
set(gcf,'renderer','painter')

%% plot incongruent responses - high ATI

%averages
prePref_Avg = arrayfun(@(x) mean(high.pre_pref{x}, 1, 'omitnan'), ...
    1:nAlign, 'UniformOutput', false);
prePref_Sem = arrayfun(@(x) sem(high.pre_pref{x}), 1:nAlign, ...
    'uniformoutput', false);

preNonpref_Avg = arrayfun(@(x) mean(high.pre_nonpref{x}, 1, 'omitnan'), ...
    1:nAlign, 'UniformOutput', false);
preNonpref_Sem = arrayfun(@(x) sem(high.pre_nonpref{x}), 1:nAlign, ...
    'uniformoutput', false);

postPref_Avg = arrayfun(@(x) mean(high.post_pref{x}, 1, 'omitnan'), ...
    1:nAlign, 'UniformOutput', false);
postPref_Sem = arrayfun(@(x) sem(high.post_pref{x}), 1:nAlign, ...
    'uniformoutput', false);

postNonpref_Avg = arrayfun(@(x) mean(high.post_nonpref{x}, 1, 'omitnan'), ...
    1:nAlign, 'UniformOutput', false);
postNonpref_Sem = arrayfun(@(x) sem(high.post_nonpref{x}), 1:nAlign, ...
    'uniformoutput', false);

%stats
pExp_post = cell(1, nAlign);
pExp_pre = pExp_post;
for jj = 1:nAlign
    pExp_post{jj} = arrayfun(@(y) permute_test(high.post_nonpref{jj}(:,y), ...
        high.post_pref{jj}(:,y), 1000, 1), xvec);
    pExp_pre{jj} = arrayfun(@(y) permute_test(high.pre_nonpref{jj}(:,y), ...
        high.pre_pref{jj}(:,y), 1000, 1), xvec);
end

% find significant bins
pExp_pre_sig = arrayfun(@(x) pExp_pre{x} < 0.05, 1:nAlign, ...
    'uniformoutput', false);
pExp_post_sig = arrayfun(@(x) pExp_post{x} < 0.05, 1:nAlign, ...
    'uniformoutput', false);

%plot
figure; hold on
tiledlayout(2, nAlign, 'TileSpacing', 'compact')

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
    scatter(x(xvec(sig)), 1, '*k')
    xline(0, '--', 'linewidth', 0.5, 'color', [0.7 0.7 0.7])
    set(gca, 'TickDir', 'out'); box off;
    axis square
    ax(ii).YRuler.TickLabelGapOffset = 1;
    title(['\rm' labels{ii}])
end
ylabel(ax(1), 'Firing rate (Hz)')
linkaxes(ax)
xlim([-0.5 1])
ylim([-0.5 1])

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
    scatter(x(xvec(sig)), 1, '*k')
    xline(0, '--', 'linewidth', 0.5, 'color', [0.7 0.7 0.7])
    set(gca, 'TickDir', 'out'); box off;
    axis square
    ax(ii).YRuler.TickLabelGapOffset = 1;
    title(['\rm' labels{ii}])
end
ylabel(ax(1), 'Firing rate (Hz)')
linkaxes(ax)
xlim([-0.5 1])
ylim([-0.4 1])

set(gcf, 'units', 'centimeters', 'position', [10 10 16 8])
set(gcf,'renderer','painter')