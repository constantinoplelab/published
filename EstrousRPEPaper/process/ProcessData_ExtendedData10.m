function ProcessData_ExtendedData10(datadir, savedir, codedir)
%ProcessData_ExtendedData10 - Process raw data saved under datadir such that it can be plotted by PlotExtendedData10.
% INPUTS:
%   datadir - Local directory where 'Golden_Nature/RawData_ExtendedData10.mat' from Zenodo was saved
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
load([datadir, 'RawData_ExtendedData10.mat'], 'OsmolalityTable', ...
    'groups')

%% Process data
%--------------------------------------------------------------------------
%ED10d. Mean serum osmolality
%--------------------------------------------------------------------------
%remove outliers
Outlier_factor = 2.5;  % outside 2.5*std is an outlier
OsmolalityTable.Osmolality(OsmolalityTable.Osmolality>mean(OsmolalityTable.Osmolality)...
    +Outlier_factor*std(OsmolalityTable.Osmolality)|OsmolalityTable.Osmolality<...
    mean(OsmolalityTable.Osmolality)-Outlier_factor*std(OsmolalityTable.Osmolality)) = NaN; %removes two samples

%%Get averages per rat
ratlist = unique(OsmolalityTable.RatID);
RatNames = cell(length(ratlist), 1);
RatGroups = cell(length(ratlist), 1);
RatConc = NaN(length(ratlist), 1);
for rat = 1:length(ratlist)
    ratT = OsmolalityTable(strcmp(OsmolalityTable.RatID, ratlist(rat)), :);
    RatNames(rat) = ratlist(rat);
    RatGroups(rat) = ratT.Group(1);
    RatConc(rat) = mean(ratT.Osmolality, 'omitnan');
end

%% Save processed data 
save([savedir 'ProcessData_ExtendedData10'], ...
    'RatConc', 'RatGroups', 'groups', '-v7.3');

end