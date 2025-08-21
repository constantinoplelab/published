function ProcessData_Figure2(datadir, savedir, codedir)
%ProcessData_Figure2 - Process raw data saved under datadir such that it can be plotted by PlotFigure2.
% INPUTS:
%   datadir - Local directory where 'Golden_Nature/RawData_Figure2.mat' from Zenodo was saved
%   codedir - Local directory where code (e.g., published/EstrousRPE/) was saved
%   savedir - Local directory where you would like the outputs to be saved

%% Set up paths
s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codedir, 'IgnoreCase', ispc);
if ~onPath % only add code dir to path if it isn't already
    addpath(genpath(codedir))
end

%% Load raw data
load([datadir, 'RawData_Figure2'],...
    'NAcc_ratlist', 'Pstruct', 'Bstruct','PstructSerum', ...
    'BstructSerum'); % Post-violation trials were excluded
load([datadir, 'RawData_Figure1'],...
    'SerumTable'); % Post-violation trials were excluded

%% Process data
%Set general variables
Stages = {'Proestrus','Diestrus'}; %estrous stages to include
Alignments = {'CPIn'}; %reward offer cue event
window = 0.5; %Measure AUC 500 ms after event
block = 1; %use mixed blocks
port = 0; %don't subset by port

%--------------------------------------------------------------------------
%2c. Heatmap of motion-corrected GRABDA signal for G024, sorted by reward volume.
%--------------------------------------------------------------------------
ratname = 'G024';
[times_resample, sorted_data_rew_rat] =...
    heatmat_byreward(Pstruct.(ratname), Bstruct.(ratname));

%--------------------------------------------------------------------------
%2d. Average response to offer cue across the 9 rats recorded from,
% separated by reward volume.
%--------------------------------------------------------------------------
[reward_da_rats, ~, T] =...
    RewardEffectAcrossRats(block, Alignments, port, NAcc_ratlist,...
    Pstruct, Bstruct);

%--------------------------------------------------------------------------
%2e.  Heatmap of motion-corrected GRABDA signal in response to reward
% offer cue for 16ul trials across over 7,000 trials recorded from example rat,
% G022, sorted by reward block.
%--------------------------------------------------------------------------
[~, sorted_data_16] = heatmat_byblock(Pstruct.(ratname),...
    Bstruct.(ratname));

%--------------------------------------------------------------------------
%2f.  Average response to 16ul offer cue across the 9 rats recorded
% from, separated by reward block.
%--------------------------------------------------------------------------
thirdrew_arg = 1;
[high_rats, low_rats] = HiloAcrossRats(thirdrew_arg, NAcc_ratlist,...
    Pstruct, Bstruct);

%--------------------------------------------------------------------------
%2g. Response to offer cue, separated by reward block and
% stage group for example rat. Gray box represents window used to
% calculate change in area under curve, 0 to 0.5 s delta AUC in h.
% Data is baseline-corrected using the 0.05 to 0 s before offer cue.
%--------------------------------------------------------------------------
ratname = 'G016';
thirdrew_arg = 0;
[high_stages_rat, low_stages_rat, err_high_stages_rat,...
    err_low_stages_rat] = hiloRew_stages(Bstruct.(ratname),...
    Pstruct.(ratname), window, thirdrew_arg, Stages);

%--------------------------------------------------------------------------
%2h. Histogram of stage effect (proestrus-diestrus) on deltaAUC for offer
% cue response, *\p<0.05, Wilcoxon signed-rank test for difference from zero,
%baseline corrected
%--------------------------------------------------------------------------
thirdrew_arg = 0;
[~, stageffect] = HiloAcrossRats_stages(window, thirdrew_arg,...
    NAcc_ratlist, Pstruct, Bstruct, Stages);

%--------------------------------------------------------------------------
%2i. Response to offer cue as a function of estradiol
%--------------------------------------------------------------------------
[delta_good, E2_good, stages] =...
    SerumEstradiol_vs_DARPE(SerumTable, PstructSerum, BstructSerum);

%--------------------------------------------------------------------------
%2j. Response to offer cue for all trials in mixed blocks, separated by reward volume and
% stage group for example rat, G037, baseline-corrected using the 0.05 to
% 0 s before offer cue.
%--------------------------------------------------------------------------
ratname = 'G016';
rewards = [1 2 3 4 5]; %standardized, linearized version of reward options
[reward_da_stages_rat, reward_da_err_stages_rat] =...
    rew_effect_stages(Bstruct.(ratname),...
    Pstruct.(ratname), block, Alignments, rewards, Stages, window);

%--------------------------------------------------------------------------
%2k. Population average response to each reward volume, separated by reward
% block, min-max normalized, baseline-corrected using the 0.05 to 0 s
% before offer cue, *p<0.05.
%--------------------------------------------------------------------------
AUC_byrat = RewardEffectByBlockAcrossRats_stages(window, Stages,...
    NAcc_ratlist, Pstruct, Bstruct);

%% Save processed data
save([savedir 'ProcessData_Figure2'],...
    'NAcc_ratlist', 'Stages', 'times_resample', 'sorted_data_rew_rat',...
    'reward_da_rats', 'T', 'sorted_data_16', 'high_rats', 'low_rats',...
    'high_stages_rat', 'low_stages_rat', 'err_high_stages_rat',...
    'err_low_stages_rat', 'stageffect', 'reward_da_stages_rat',...
    'reward_da_err_stages_rat', 'AUC_byrat', ...
    'delta_good', 'E2_good', 'stages');

end