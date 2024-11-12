function plotFigure4(ephysPath_expert, incongruentPath_expert, incongruentPath_naive)
% plot data from figure 4
% YOU MUST RUN THE FUNCTION processIncongruentResponseData FOR BOTH EXPERT
% AND NAIVE DATA TO RUN THIS FUNCTION

%pref and nonpref nomenclature is based on preferred and non-preferred
%adaptation block (opposite of transition type naming in the paper)

%INPUTS:
%   ephysPath_expert = path to ephys data structs for expert rats (e.g.
%       'E003_2022-08-01.mat')
%   incongruentPath_expert = path to file saved after running
%       processIncongruentResponseData on expert ephys data
%   incongruentPath_naive = path to files saved after running
%       processIncongruentResponseData on expert ephys data


%% general

alignto = {'COFF' 'SON' 'SOFF' 'Rew'};
nAlign = length(alignto);

timeBin = [-4 8]; % time window around each event saved in the SU structs
x = timeBin(1):0.05:timeBin(2);
wndw = [0 0.5]; %window for determining significance
w = arrayfun(@(y) find(x == wndw(y)), 1:2); %indices of window for computing significance

sem = @(x) std(x,'omitnan')./sqrt(sum(any(~isnan(x), 2)));

setDefaultFigProps

%% process incongruent responses

%incongruent responses for experts
expert = load([incongruentPath_expert, 'incongruentResponses.mat']);

%incongruent responses for naive
naive = load([incongruentPath_naive, 'incongruentResponses.mat']);

%% example block significant cell
filename = 'S021_2021-09-22.mat';
cellID = 150;

load(strcat(ephysPath_expert, filename), 'S', 'SU');
ids = extractfield(SU, 'cluster_id');
use = ids == cellID;

hmat = SU(use).hmat.Rew;
byBlock = arrayfun(@(x) mean(hmat(S.Block == x, :), 1, 'omitnan'), ...
    [2, 3], 'UniformOutput', false);
byBlock_sem = arrayfun(@(x) sem(hmat(S.Block == x, :)), [2 3], ...
    'UniformOutput', false);

%% average incongruent and congruent trial responses

%Experts
incongPref_Avg = arrayfun(@(x) mean(expert.incongPref{x}, 1, 'omitnan'), ...
    1:nAlign, 'UniformOutput', false);
incongPref_Sem = arrayfun(@(x) sem(expert.incongPref{x}), 1:nAlign, ...
    'uniformoutput', false);

incongNonpref_Avg = arrayfun(@(x) mean(expert.incongNonpref{x}, 1, 'omitnan'), ...
    1:nAlign, 'UniformOutput', false);
incongNonpref_Sem = arrayfun(@(x) sem(expert.incongNonpref{x}), 1:nAlign, ...
    'uniformoutput', false);

congPref_Avg = arrayfun(@(x) mean(expert.congPref{x}, 1, 'omitnan'), ...
    1:nAlign, 'UniformOutput', false);
congPref_Sem = arrayfun(@(x) sem(expert.congPref{x}), 1:nAlign, ...
    'uniformoutput', false);

congNonpref_Avg = arrayfun(@(x) mean(expert.congNonpref{x}, 1, 'omitnan'), ...
    1:nAlign, 'UniformOutput', false);
congNonpref_Sem = arrayfun(@(x) sem(expert.congNonpref{x}), 1:nAlign, ...
    'uniformoutput', false);


%Naive
incongPref_Avg_N = arrayfun(@(x) mean(naive.incongPref{x}, 1, 'omitnan'), ...
    1:nAlign, 'UniformOutput', false);
incongPref_Sem_N = arrayfun(@(x) sem(naive.incongPref{x}), 1:nAlign, ...
    'uniformoutput', false);

incongNonpref_Avg_N = arrayfun(@(x) mean(naive.incongNonpref{x}, 1, 'omitnan'), ...
    1:nAlign, 'UniformOutput', false);
incongNonpref_Sem_N = arrayfun(@(x) sem(naive.incongNonpref{x}), 1:nAlign, ...
    'uniformoutput', false);

congPref_Avg_N = arrayfun(@(x) mean(naive.congPref{x}, 1, 'omitnan'), ...
    1:nAlign, 'UniformOutput', false);
congPref_Sem_N = arrayfun(@(x) sem(naive.congPref{x}), 1:nAlign, ...
    'uniformoutput', false);

congNonpref_Avg_N = arrayfun(@(x) mean(naive.congNonpref{x}, 1, 'omitnan'), ...
    1:nAlign, 'UniformOutput', false);
congNonpref_Sem_N = arrayfun(@(x) sem(naive.congNonpref{x}), 1:nAlign, ...
    'uniformoutput', false);

%% compute p-values using a non-parametric permutation test
%reported p-values are bonferroni corrected for 4 comparisons. values may 
%vary slightly based on the permutation
 
%stats
bonferroniCorrect = 4;

inNp_E = arrayfun(@(x) mean(expert.incongNonpref{x}(:, w(1):w(2)), 2, 'omitnan'), ...
    1:nAlign, 'uniformoutput', false);
inNp_E = arrayfun(@(x) inNp_E{x}(~isnan(inNp_E{x})), 1:nAlign, ...
    'uniformoutput', false);
inNp_N = arrayfun(@(x) mean(naive.incongNonpref{x}(:, w(1):w(2)), 2, 'omitnan'), ...
    1:nAlign, 'uniformoutput', false);
inNp_N = arrayfun(@(x) inNp_N{x}(~isnan(inNp_N{x})), 1:nAlign, ...
    'uniformoutput', false);
[p_inNp, ~] = arrayfun(@(x) permute_test(inNp_N{x}, inNp_E{x}, 1000), ...
    1:nAlign, 'uniformoutput', false); 

inP_E = arrayfun(@(x) mean(expert.incongPref{x}(:, w(1):w(2)), 2, 'omitnan'), ...
    1:nAlign, 'uniformoutput', false);
inP_E = arrayfun(@(x) inP_E{x}(~isnan(inP_E{x})), 1:nAlign, ...
    'uniformoutput', false);
inP_N = arrayfun(@(x) mean(naive.incongPref{x}(:, w(1):w(2)), 2, 'omitnan'), ...
    1:nAlign, 'uniformoutput', false);
inP_N = arrayfun(@(x) inP_N{x}(~isnan(inP_N{x})), 1:nAlign, ...
    'uniformoutput', false);
[p_inP, ~] = arrayfun(@(x) permute_test(inP_N{x}, inP_E{x}, 1000), ...
    1:nAlign, 'uniformoutput', false); 

cNp_E = arrayfun(@(x) mean(expert.congNonpref{x}(:, w(1):w(2)), 2, 'omitnan'), ...
    1:nAlign, 'uniformoutput', false);
cNp_E = arrayfun(@(x) cNp_E{x}(~isnan(cNp_E{x})), 1:nAlign, ...
    'uniformoutput', false);
cNp_N = arrayfun(@(x) mean(naive.congNonpref{x}(:, w(1):w(2)), 2, 'omitnan'), ...
    1:nAlign, 'uniformoutput', false);
cNp_N = arrayfun(@(x) cNp_N{x}(~isnan(cNp_N{x})), 1:nAlign, ...
    'uniformoutput', false);
[p_cNp, ~] = arrayfun(@(x) permute_test(cNp_N{x}, cNp_E{x}, 1000), ...
    1:nAlign, 'uniformoutput', false); 

cP_E = arrayfun(@(x) mean(expert.congPref{x}(:, w(1):w(2)), 2, 'omitnan'), ...
    1:nAlign, 'uniformoutput', false);
cP_E = arrayfun(@(x) cP_E{x}(~isnan(cP_E{x})), 1:nAlign, ...
    'uniformoutput', false);
cP_N = arrayfun(@(x) mean(naive.congPref{x}(:, w(1):w(2)), 2, 'omitnan'), ...
    1:nAlign, 'uniformoutput', false);
cP_N = arrayfun(@(x) cP_N{x}(~isnan(cP_N{x})), 1:nAlign, ...
    'uniformoutput', false);
[p_cP, ~] = arrayfun(@(x) permute_test(cP_N{x}, cP_E{x}, 1000), ...
    1:nAlign, 'uniformoutput', false); 

%% plot

figure; hold on
tiledlayout(2, 5, 'TileSpacing', 'tight')

ax1 = nexttile;
shadedErrorBar(x, byBlock{1}, byBlock_sem{1}, 'lineprops', {'color', ...
    'r', 'linewidth', 1})
shadedErrorBar(x, byBlock{2}, byBlock_sem{2}, 'lineprops', {'color', ...
    'b', 'linewidth', 1})
set(gca, 'TickDir', 'out'); box off;
xlim([-0.5 2]);
xlabel('Time from reward (s)')
ylabel('Firing rate (Hz)')
axis square
ax1.YRuler.TickLabelGapOffset = 1;
title('\rm Block cell')

ax2 = nexttile;
shadedErrorBar(x, incongNonpref_Avg{4}, incongNonpref_Sem{4}, 'lineprops', ...
    {'color', 'k', 'linewidth', 1})
shadedErrorBar(x, incongNonpref_Avg_N{4}, incongNonpref_Sem_N{4}, 'lineprops', ...
    {'color', 'm', 'linewidth', 1})
text(0, 4.5, string(p_inNp{4}*bonferroniCorrect))
set(gca, 'TickDir', 'out'); box off;
axis square
ax2.YRuler.TickLabelGapOffset = 1;
title('\rm Incongruent nonpref')

ax3 = nexttile;
shadedErrorBar(x, incongPref_Avg{4}, incongPref_Sem{4}, 'lineprops', ...
    {'color', 'k', 'linewidth', 1})
shadedErrorBar(x, incongPref_Avg_N{4}, incongPref_Sem_N{4}, 'lineprops', ...
    {'color', 'm', 'linewidth', 1})
text(0, 4.5, string(p_inP{4}*bonferroniCorrect))
set(gca, 'TickDir', 'out'); box off;
axis square
ax3.YRuler.TickLabelGapOffset = 1;
title('\rm Incongruent')

ax4 = nexttile;
shadedErrorBar(x, congNonpref_Avg{4}, congNonpref_Sem{4}, 'lineprops', ...
    {'color', 'k', 'linewidth', 1})
shadedErrorBar(x, congNonpref_Avg_N{4}, congNonpref_Sem_N{4}, 'lineprops', ...
    {'color', 'm', 'linewidth', 1})
text(0, 4.5, string(p_cNp{4}*bonferroniCorrect))
set(gca, 'TickDir', 'out'); box off;
axis square
ax4.YRuler.TickLabelGapOffset = 1;
title('\rm Congruent nonpref')

ax5 = nexttile;
shadedErrorBar(x, congPref_Avg{4}, congPref_Sem{4}, 'lineprops', ...
    {'color', 'k', 'linewidth', 1})
shadedErrorBar(x, congPref_Avg_N{4}, congPref_Sem_N{4}, 'lineprops', ...
    {'color', 'm', 'linewidth', 1})
text(0, 4.5, string(p_cP{4}*bonferroniCorrect))
set(gca, 'TickDir', 'out'); box off;
axis square
ax5.YRuler.TickLabelGapOffset = 1;
title('\rm Congruent')

%blank spot for schematic
nexttile
axis square

labels = {'Trial start' 'LED on' 'LED off' 'Reward'};
for ii = 1:nAlign

    ax(ii) = nexttile;
    shadedErrorBar(x, incongNonpref_Avg{ii}, incongNonpref_Sem{ii}, 'lineprops', ...
        {'color', 'k', 'linewidth', 1})
    shadedErrorBar(x, incongNonpref_Avg_N{ii}, incongNonpref_Sem_N{ii}, 'lineprops', ...
        {'color', 'm', 'linewidth', 1})
    text(0, 4.5, string(p_inNp{ii}*bonferroniCorrect))
    set(gca, 'TickDir', 'out'); box off;
    axis square
    ax(ii).YRuler.TickLabelGapOffset = 1;
    title(['\rm' labels{ii}])
end

linkaxes([ax2, ax3, ax4, ax5 ax])
xlim([-1 3])
ylim([4 14])
ylabel(ax2, 'Firing rate (Hz)')
xlabel(ax4, 'Time from reward (s)')
set(gcf, 'units', 'centimeters', 'position', [10 10 20 8])
set(gcf,'renderer','painter')


