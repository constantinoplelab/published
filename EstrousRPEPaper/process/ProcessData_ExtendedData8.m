function ProcessData_ExtendedData8(datadir, savedir, codedir)
%ProcessData_ExtendedData8 - Process raw data saved under datadir such that it can be plotted by PlotExtendedData8.
% INPUTS:
%   datadir - Local directory where 'Golden_Nature/RawData_ExtendedData8.mat' from Zenodo was saved
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
load([datadir, 'RawData_ExtendedData8.mat'], 'results', 'ratlist')

%% Process data
%--------------------------------------------------------------------------
% ED8c. TH and ChR2 overlap
%--------------------------------------------------------------------
percentages = [];
for f=1:length(ratlist)
    ratname = ratlist{f};
    percentages = [percentages; results.(ratname){:,4}];
end

%% Save processed data 
save([savedir 'ProcessData_ExtendedData8'], 'percentages', '-v7.3');

end