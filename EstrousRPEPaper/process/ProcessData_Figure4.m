function ProcessData_Figure4(datadir, codedir, savedir)
%ProcessData_Figure4 - Process raw data saved under datadir such that it can be plotted by PlotFigure4.
% INPUTS:
%   datadir - Local directory where 'Golden_Nature/RawData_Figure4.mat' from Zenodo was saved
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
load([datadir, 'RawData_Figure4'],...
    'Pstruct', 'Bstruct', 'ratTrial', 'Aligned_rat', 'S',...
    'optoratlist', 'isoptorat', 'optoafterdates', 'AlignedRats', 'SRats')

%% Process data
%Set general variables
NAcc_ratlist = {'G008','G016','G021','G022','G024',...
    'G036','G037','G051','G127'};

%--------------------------------------------------------------------------
%4a-c. DA trace at offer cue by following initiation time bin. No
%baseline-correction, mixed block and violation trials only.
%--------------------------------------------------------------------------
%initiation time histogram inset
binnum = 7;
[photITIs, ITIbins] =...
    get_ITIbins(NAcc_ratlist, Bstruct, binnum);
%DA signal
event = 'CPIn';
window = 0.5;
[da_rats, AUC_overrats, P_byrat, Rcorr_byrat, bins, T, ITI_legend] =...
    ITIEffectAcrossRats(NAcc_ratlist, ITIbins, binnum, event, window,...
    Pstruct, Bstruct); 

%--------------------------------------------------------------------------
%4g. Model of initiation time with reinforcement learning model,
% simulation of opto stimulation with additional RPE on 30% of trials
%--------------------------------------------------------------------------
nexttile
opto_args = [0 1];
LatOpts = cell(2, 1);
for oa=1:2
    [LatOpts{oa}, ~, ~, ~] =...
        GenerateLatencyData_VanillaAlpha_cg(0.35, ratTrial,...
        opto_args(oa), 0, 4);
end

%--------------------------------------------------------------------------
%4h. Example opto rat, G073
%--------------------------------------------------------------------------
opto_start_date = '20230214';
[avgITIbyblock_opto, errITIbyblock_opto,...
    avgITIbyblock_nonopto, errITIbyblock_nonopto] =...
    CompareOptoITIbyBlock(opto_start_date, Aligned_rat, S);

%--------------------------------------------------------------------------
%4i. Average effect of opto on initiation time for TH-Cre + ChR2 rats
% (opto compared to control sessions)
%--------------------------------------------------------------------------
[nonopto_norm, opto_norm] =...
    CompareOptoITIbyBlock_overrats(optoratlist, optoafterdates,...
    AlignedRats, SRats);

%% Save processed data
save([savedir 'ProcessData_Figure4'], 'NAcc_ratlist', 'photITIs', 'ITIbins',...
    'da_rats', 'AUC_overrats', 'P_byrat', 'Rcorr_byrat', 'bins', 'T',...
    'ITI_legend', 'LatOpts', 'avgITIbyblock_opto', 'errITIbyblock_opto',...
    'avgITIbyblock_nonopto', 'errITIbyblock_nonopto', 'nonopto_norm',...
    'opto_norm', 'isoptorat', 'ratTrial');

end