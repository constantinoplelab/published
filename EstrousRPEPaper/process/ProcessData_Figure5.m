function ProcessData_Figure5(datadir, codedir, savedir)
%ProcessData_Figure5 - Process raw data saved under datadir such that it can be plotted by PlotFigure5.
% INPUTS:
%   datadir - Local directory where 'Golden_Nature/RawData_Figure5.mat' from Zenodo was saved
%   codedir - Local directory where code (e.g., published/EstrousRPE/) was saved
%   savedir - Local directory where you would like the outputs to be saved

%% Set up paths
s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codedir, 'IgnoreCase', ispc);
if ~onPath % only add code dir to path if it already isn't
    addpath(genpath(codedir))
end

%% Load raw data
load([datadir, 'RawData_Figure5'],...
    'MS_pvd', 'MS_evd', 'TABLE_SERT', 'TABLE_DAT', 'EMdata')

%% Process data

%--------------------------------------------------------------------------
%5a. Volcano plots
%--------------------------------------------------------------------------
comparison = 1;
direction = -1;
%proestrus/diestrus
[log2FC_pvd, y_pvd, txt1_pvd, txt2_pvd, txt3_pvd,...
    stagescompared_pvd, stageslfc_pvd] = volcanoplot(MS_pvd, comparison,...
    direction);

%estrus/diestrus
comparison = 3;
direction = -1;
[log2FC_evd, y_evd, txt1_evd, txt2_evd, txt3_evd,...
    stagescompared_evd, stageslfc_evd] = volcanoplot(MS_evd, comparison,...
    direction);

%--------------------------------------------------------------------------
%5c. SERT pixel area
%--------------------------------------------------------------------------
[Stages_SERT, AreaByStage_SERT, ~, pvals_SERT] =...
    SERTDATquantplot(TABLE_SERT, 'SERT');

%--------------------------------------------------------------------------
%5d. DAT pixel area
%--------------------------------------------------------------------------
[Stages_DAT, AreaByStage_DAT, ~, pvals_DAT] =...
    SERTDATquantplot(TABLE_DAT, 'DAT');

%--------------------------------------------------------------------------
%5f. Membranous vs intracellular DAT particles identified with electron
%microscopy and their distance from the membrane
%--------------------------------------------------------------------------
[EMdata, FracSIGIC, FracSIGMem, groups, groupID] =...
    DATlocation(EMdata);


%% Save processed data
save([savedir 'ProcessData_Figure5'], 'MS_pvd', 'log2FC_pvd', 'y_pvd',...
    'txt1_pvd', 'txt2_pvd', 'txt3_pvd', 'stagescompared_pvd',...
    'stageslfc_pvd', 'MS_evd', 'log2FC_evd', 'y_evd',...
    'txt1_evd', 'txt2_evd', 'txt3_evd', 'stagescompared_evd',...
    'stageslfc_evd', 'Stages_SERT', 'AreaByStage_SERT', 'pvals_SERT',...
    'Stages_DAT', 'AreaByStage_DAT', 'pvals_DAT', 'EMdata',...
    'FracSIGIC', 'FracSIGMem', 'groups', 'groupID');

end