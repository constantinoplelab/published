function ProcessData_ExtendedData3(datadir, savedir, codedir)
%ProcessData_ExtendedData3 - Process raw data saved under datadir such that it can be plotted by PlotExtendedData3.
% INPUTS:
%   datadir - Local directory where 'Golden_Nature/RawData_Figure1.mat' from Zenodo was saved
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
%Rat lists
load([datadir, 'RawData_Figure1'], 'RatBehaviorData');

%% Create variables
NAcc_ratlist = {'G008','G016','G021','G022','G024',...
    'G036','G037','G051','G127','G138','A007','A008','L060'}; 
Serum_ratlist = {'G125','G127','G128','G131','G133','G134','G135',...
    'G136','G138','G140','G141','A007','A008','L060','Q001','Q002',...
    'Q003'}; 

%% Process data
%Initiation times by block for specific rat groups
% Get data over estrous
ITIbyBlockPhotometry = ITIbyBlock_estrous(NAcc_ratlist, RatBehaviorData); 
ITIbyBlockPhotometry_fem = ITIbyBlock_females(NAcc_ratlist, RatBehaviorData);
ITIbyBlockPhotometry.Female = ITIbyBlockPhotometry_fem.Female;

ITIbyBlockSerum = ITIbyBlock_estrous(Serum_ratlist, RatBehaviorData);
ITIbyBlockSerum_fem = ITIbyBlock_females(Serum_ratlist, RatBehaviorData);
ITIbyBlockSerum.Female = ITIbyBlockSerum_fem.Female;

%% Save processed data
save([savedir 'ProcessData_ExtendedData3'], 'NAcc_ratlist', 'Serum_ratlist',...
    'ITIbyBlockPhotometry', 'ITIbyBlockSerum')

end