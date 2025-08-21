function ProcessData_ExtendedData2(datadir, savedir, codedir)
%ProcessData_ExtendedData2 - Process raw data saved under datadir such that it can be plotted by PlotExtendedData2.
% INPUTS:
%   datadir - Local directory where 'Golden_Nature/RawData_ExtendedData2.mat' from Zenodo was saved
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
load([datadir, 'RawData_Figure1'], 'RatBehaviorData', 'f_ratlist');

%% Process data
WTbyBlock = WTbyBlock_stages(f_ratlist, RatBehaviorData);

%% Save processed data
save([savedir 'ProcessData_ExtendedData2'], 'WTbyBlock', 'f_ratlist');

end