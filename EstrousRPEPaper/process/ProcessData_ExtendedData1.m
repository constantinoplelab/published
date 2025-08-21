function ProcessData_ExtendedData1(datadir, savedir, codedir)
%ProcessData_ExtendedData1 - Process raw data saved under datadir such that it can be plotted by PlotExtendedData1.
% INPUTS:
%   datadir - Local directory where 'Golden_Nature/RawData_ExtendedData1.mat' from Zenodo was saved
%   codedir - Local directory where code (e.g., published/EstrousRPE/) was saved
%   savedir - Local directory where you would like the outputs to be saved
%   and where ProcessData_Figure1 was saved

%% Set up paths
s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codedir, 'IgnoreCase', ispc);
if ~onPath % only add code dir to path if it isn't already
    addpath(genpath(codedir))
end

%% Load raw data
load([savedir, 'ProcessData_Figure1'], 'f_ratlist',...
    'ITIbyBlock', 'beh_sens_staged');
load([datadir, 'RawData_Figure1'], 'SerumTable');
load([datadir, 'RawData_ExtendedData1'], 'TrunkBloodEstradiolTable');

%% Process data
%Set general variables
Stages = {'Metestrus', 'Diestrus', 'Proestrus', 'Estrus'};

%--------------------------------------------------------------------------
%ED1c. Estradiol serum levels, with rats as individual points
%--------------------------------------------------------------------------
proestrus = TrunkBloodEstradiolTable(strcmp(TrunkBloodEstradiolTable.Stage, 'proestrus'), :); %subset samples predicted to be in proestrus
diestrus = TrunkBloodEstradiolTable(strcmp(TrunkBloodEstradiolTable.Stage, 'diestrus'), :); %subset samples predicted to be in diestrus
estrus = TrunkBloodEstradiolTable(strcmp(TrunkBloodEstradiolTable.Stage, 'estrus'), :); %subset samples predicted to be in estrus
%Concentration split by stage
%split into early/late proestrus
earlyproidx = proestrus.Time < 0.6;
lateproidx = proestrus.Time > 0.6;
Concentration.D = diestrus.Concentration;
Concentration.PE = proestrus.Concentration(earlyproidx);
Concentration.PL = proestrus.Concentration(lateproidx);
Concentration.E = estrus.Concentration;


%% Save processed data
save([savedir 'ProcessData_ExtendedData1'], 'f_ratlist',...
    'Stages', 'Concentration', 'SerumTable','ITIbyBlock',...
    'beh_sens_staged');

end