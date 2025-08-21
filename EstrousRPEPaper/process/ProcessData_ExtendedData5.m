function ProcessData_ExtendedData5(datadir, savedir, codedir)
%ProcessData_ExtendedData5 - Process raw data saved under datadir such that it can be plotted by PlotExtendedData5.
% INPUTS:
%   datadir - Local directory where 'Golden_Nature/RawData_ExtendedData5.mat' from Zenodo was saved
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
    'NAcc_ratlist'); 
load([datadir, 'RawData_Figure3'], 'Bstruct', 'PstructCPIn',...
    'Proestrus_last', 'Estrus_last', 'Metestrus_last', 'Diestrus_last',...
    'pro_ratList', 'est_ratList', 'met_ratList', 'di_ratList') % Post-violation trials were excluded
load([datadir, 'RawData_ExtendedData5'], 'mCherryG016', 'GRABDAG016',...
    'bstructG016')

%% Process data
%Set general variables
Stages = {'Metestrus','Diestrus','Proestrus','Estrus'}; %estrous stages to include
window = 0.5; %Measure AUC 500 ms after event

%--------------------------------------------------------------------------
%ED5a-b. Min-max normalized AUC for the largest reward (64 ul) during 
% mixed and high blocks over the cycle
%--------------------------------------------------------------------------
disp('getting AUC by reward volume')
AUC_byrat = RewardEffectByBlockAcrossRats_stages(window,...
    Stages, NAcc_ratlist, PstructCPIn, Bstruct);

%--------------------------------------------------------------------------
%ED5c. Min-max normalized AUC for the largest RPE over the cycle
%--------------------------------------------------------------------------
numbins = 6;
window = 0.5; %AUC window
event = 'CPIn'; %align to reward offer cue
disp('getting AUC by RPE')
[pro_DA_binned, est_DA_binned, met_DA_binned, di_DA_binned,...
    RPEbins] = DA_by_RPE_estrous(NAcc_ratlist,...
    PstructCPIn, Bstruct, Proestrus_last, Estrus_last, Metestrus_last,...
    Diestrus_last, pro_ratList, est_ratList, met_ratList, di_ratList,...
    numbins, window, event); %all blocks, last 20 trials, equally spaced bins

%--------------------------------------------------------------------------
%ED5d-e. Motion-corrected GRABDA and mCherry split by reward during mixed 
% blocks, aligned to each event
%--------------------------------------------------------------------------
ratname = 'G016';
port = 0;
rewards = [1 2 3 4 5];
Alignments = {'CPOn', 'CPIn', 'SideOn', 'SideOff', 'Reward', 'OptOut'};
[RewardDA_gfp, RewardDA_err_gfp, T, ~] =...
    rew_effect(bstructG016.bstruct.(ratname),...
    GRABDAG016.pstruct.(ratname), 1, rewards, Alignments,...
    port);
[RewardDA_mCherry, RewardDA_err_mCherry, ~, ~] =...
    rew_effect(bstructG016.bstruct.(ratname),...
    mCherryG016.pstruct.(ratname), 1, rewards, Alignments,...
    port); 

%--------------------------------------------------------------------------
%ED5f-g. Motion-corrected GRABDA and mCherry split by block, aligned to each event
%--------------------------------------------------------------------------
[hi_gfp, lo_gfp, high_err_gfp, low_err_gfp] =...
    hilo_all_alignments(ratname, bstructG016.bstruct,...
    GRABDAG016.pstruct.(ratname), 0, 0);
[hi_mCherry, lo_mCherry, high_err_mCherry, low_err_mCherry] =...
    hilo_all_alignments(ratname, bstructG016.bstruct,...
    mCherryG016.pstruct.(ratname), 0, 0);


%% Save processed data
save([savedir 'ProcessData_ExtendedData5'],...
    'NAcc_ratlist', 'Stages', 'AUC_byrat', 'pro_DA_binned',...
    'est_DA_binned', 'met_DA_binned', 'di_DA_binned', 'RPEbins',...
    'RewardDA_gfp', 'RewardDA_err_gfp', 'T', 'RewardDA_mCherry',...
    'RewardDA_err_mCherry', 'hi_gfp', 'lo_gfp', 'high_err_gfp',...
    'low_err_gfp', 'hi_mCherry', 'lo_mCherry', 'high_err_mCherry',...
    'low_err_mCherry', 'Alignments', 'rewards');

end