function ProcessData_ExtendedData7(datadir, savedir, codedir)
%ProcessData_ExtendedData7 - Process raw data saved under datadir such that it can be plotted by PlotExtendedData7.
% INPUTS:
%   datadir - Local directory where 'Golden_Nature/RawData_ExtendedData7.mat' from Zenodo was saved
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
load([datadir, 'RawData_ExtendedData7.mat'], 'NAcc_ratlist', ...
    'Bstruct', 'PstructSOn', 'PstructSOff')

%% Process data
%--------------------------------------------------------------------------
% ED7a. Dopamine during the delay, split by delay length.
% Includes negative ramping, which suggests state inference.
%--------------------------------------------------------------------
numbins = 6;
delaybins = get_delays(NAcc_ratlist, Bstruct, numbins);    
event = 'SideOn';
block = 1; %Mixed block only
[da_rats_SideOn, T, bins, delay_legend] = DelayEffectAcrossRats(NAcc_ratlist,...
    block, delaybins, event, Bstruct, PstructSOn);

%--------------------------------------------------------------------------
% ED7b. Dopamine at the time of the side light turning off (cues that the
% delay is over and the reward is available), split by delay length.NAcc_ratlist
%--------------------------------------------------------------------------
event = 'SideOff';
da_rats_SideOff = DelayEffectAcrossRats(NAcc_ratlist,...
    block, delaybins, event, Bstruct, PstructSOff);

%--------------------------------------------------------------------------
% ED7c-d. Phasic response at the time of the side light turning off, split by
% estrous stage (additional evidence that RPE is enhanced in proestrus).
% Time course and quantification of area under the curve.
%--------------------------------------------------------------------------
window = 0.5; %measure area under the curve for the first 500 ms after the event
Stages = {'Proestrus','Diestrus'};      
[da_rats_stages, AUC_stages] =...
    DelayEffectAcrossRats_stages(NAcc_ratlist, block, delaybins,...
    window, PstructSOff, Bstruct, event, Stages);

%% Save processed data 
save([savedir 'ProcessData_ExtendedData7'],'NAcc_ratlist','numbins', ...
    'da_rats_SideOn','T', 'bins', 'delay_legend', 'da_rats_SideOff', ...
    'Stages', 'da_rats_stages', 'AUC_stages');

end