function ProcessData_ExtendedData4(datadir, savedir, codedir)
%ProcessData_ExtendedData4 - Process raw data saved under datadir such that it can be plotted by PlotExtendedData4.
% INPUTS:
%   datadir - Local directory where 'Golden_Nature/RawData_ExtendedData4.mat' from Zenodo was saved
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

%% Create variables
NAcc_ratlist = {'G008','G016','G021','G022','G024',...
    'G036','G037','G051','G127','G138','A007','A008','L060'}; 
rewards = 1:5;
Alignments = {'CPOn', 'CPIn', 'SideOn', 'SideOff', 'Reward', 'OptOut'};
port = 0;

%% Process data
DAbyrat = cell(length(NAcc_ratlist), 1);
DAerrbyrat = cell(length(NAcc_ratlist), 1);
for rat = 1:length(NAcc_ratlist)
    ratname = NAcc_ratlist{rat};
    disp(ratname)

    %Load rat data
    load([datadir 'RawData_' ratname '_ExtendedData4'], 'bstruct', 'pstruct');

    [DAbyrat{rat}, DAerrbyrat{rat}, T] = rew_effect(bstruct,...
        pstruct, 1, rewards, Alignments, port); %rmv_postvios
end

%% Save processed data
save([savedir 'ProcessData_ExtendedData4'], 'NAcc_ratlist', 'rewards',...
    'Alignments', 'DAbyrat', 'DAerrbyrat', 'T')

end