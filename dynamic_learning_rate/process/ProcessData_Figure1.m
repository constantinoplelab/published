function ProcessData_Figure1(BestFitEarly, codedir, savedir)
%ProcessData_Figure1 - function to process data to generate Figure 1. Saves
% data in savedir to be accessed by PlotFigure1
% INPUTS:
% BestFitEarly - BestFit.mat file found at: Mah_CellReports_BehavorialData/FitAll_ITI_VanillaAlpha_First10/BestFit.mat
%   from Zenodo
% codedir - Directory of code (e.g., published/dynamic_learning_rate/)
% savedir - Location you would like the outputs to be saved

%% Set up paths

s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codedir, 'IgnoreCase', ispc);
if ~onPath % only add code dir to path if it already isn't
    addpath(genpath(codedir))
end

%% Process Data

% Pull list of rats 
fieldlist = fields(BestFitEarly);
ratList = fieldlist(structfun(@isstruct, BestFitEarly));

nback = 10; % Number of previous trials to regress over

exampleRat = 'J031'; % Example individual rat identifier
exampleRatI = strcmp(ratList, exampleRat);

% Preallocate data structures
itiByBlk = struct('raw', nan(length(ratList), 6),...
    'z', nan(length(ratList), 6)); % ITI as a function of block
pITIbyBlkRat = nan(length(ratList), 3); % P-values

betasAll = nan(length(ratList), nback+1); % Previous reward regression

% Loop over rats
for rr = 1:length(ratList)
    fprintf('%d out of %d\n', rr, length(ratList)) % Check progress

    A = BestFitEarly.(ratList{rr}).All.ratTrial; % Pull rat's data

    % ITI by block - rat
    [hi, lo, mi, pITIbyBlkRat(rr,:)] = iticurves(A);
    itiByBlk.raw(rr,:) = [lo.raw, mi.raw, hi.raw];
    itiByBlk.z(rr,:) = [lo.z, mi.z, hi.z];

    % Previous reward regression
    betasAll(rr,:) =...
        regress_latency_vs_rew(A, nback, false, true, true, true);
end

% P-value for population ITI as a function of block 
pITIbyBlkPop = [signrank(itiByBlk.z(:,1), itiByBlk.z(:,3));
    signrank(itiByBlk.z(:,1), itiByBlk.z(:,2));
    signrank(itiByBlk.z(:,2), itiByBlk.z(:,3))];

% P-value for previous regression coefficients
pRegressionPop =...
    arrayfun(@(i) signrank(betasAll(:,i)), 1:nback+1);

% Save Data
if nargin == 3
    save([savedir, 'Figure1_Data'], 'ratList', 'exampleRat',...
        'exampleRatI', 'itiByBlk', 'pITIbyBlkRat', 'betasAll',...
        'pITIbyBlkPop', 'pRegressionPop')
end

end