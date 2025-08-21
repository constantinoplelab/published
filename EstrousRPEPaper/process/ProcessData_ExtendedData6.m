function ProcessData_ExtendedData6(datadir, savedir, codedir)
%ProcessData_ExtendedData6 - Process raw data saved under datadir such that it can be plotted by PlotExtendedData6.
% INPUTS:
%   datadir - Local directory where 'Golden_Nature/RawData_ExtendedData6.mat' from Zenodo was saved
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
load([datadir, 'RawData_ExtendedData6.mat'], 'CPu_ratlist', 'Pstruct_CPu', ...
    'Bstruct_CPu')

%% Process data
%Set general variables
window = 0.5;
Stages = {'Proestrus', 'Diestrus'};
rewards = [4 8 16 32 64];

%--------------------------------------------------------------------------
% ED6d. Event aligned DA response, split by reward in mixed block for off-target rats.
%--------------------------------------------------------------------------
block = 1; %mixed block
port = 0;
Alignments = {'CPOn', 'CPIn', 'SideOn', 'SideOff', 'Reward', 'OptOut'};
[da_rats_cpu, ~, T] = RewardEffectAcrossRats(block, Alignments, port, ...
    CPu_ratlist, Pstruct_CPu, Bstruct_CPu); %remove post violations

%--------------------------------------------------------------------------
% ED6e. Response to offer cue for all trials in mixed blocks, separated by reward volume and
% stage group for example rat, G027, baseline-corrected using the 0.05 to
% 0 s before offer cue.
%--------------------------------------------------------------------------
ratname = 'G027';
[DA_G027, DA_err_G027] = rew_effect_stages(Bstruct_CPu.(ratname),...
    Pstruct_CPu.(ratname), block, Alignments(2), rewards, Stages, window);

%--------------------------------------------------------------------------
% ED6f. Response to offer cue, separated by reward block and
% stage group for example rat, G027. Gray box represents window used to
% calculate change in area under curve, 0 to 0.5 s delta AUC in g.
% Data is baseline-corrected using the 0.05 to 0 s before offer cue.
%--------------------------------------------------------------------------
thirdrew_arg = 0; %pool all rewards in each block
[high_cpu, low_cpu, err_high, err_low] = hiloRew_stages(Bstruct_CPu.(ratname),...
    Pstruct_CPu.(ratname), window, thirdrew_arg, Stages);

%--------------------------------------------------------------------------
% ED6g. Population average response to each reward volume, separated by reward
% block, min-max normalized, baseline-corrected using the 0.05 to 0 s
% before offer cue, *p<0.05.
%--------------------------------------------------------------------------
AUC_byrat_CPu =...
    RewardEffectByBlockAcrossRats_stages(window, Stages, CPu_ratlist,...
    Pstruct_CPu, Bstruct_CPu);

%--------------------------------------------------------------------------
% ED6h. Histogram of stage effect (fertile-non-fertile) on deltaAUC for offer
% cue response to block (high - low)
%--------------------------------------------------------------------------
thirdrew_arg = 0; %pool all rewards in each block
[~, stageffect_CPu] = HiloAcrossRats_stages(window, thirdrew_arg,...
    CPu_ratlist, Pstruct_CPu, Bstruct_CPu, Stages);

%% Save processed data
save([savedir 'ProcessData_ExtendedData6'],'Alignments', 'rewards',...
    'da_rats_cpu','T','Stages','DA_G027','DA_err_G027',...
    'high_cpu', 'low_cpu', 'err_high', 'err_low', 'AUC_byrat_CPu',...
    'stageffect_CPu', 'CPu_ratlist');

end