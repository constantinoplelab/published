function FigureS10(datadir, codedir)
%FigureS10 -  Sub-optimal inferential model with lambda
%   datadir = directory of dataset
%   codedir = director of code

%% Set up paths
s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codedir, 'IgnoreCase', ispc);

if ~onPath
    addpath(genpath(codedir))
end

%% Load data

mdl3 = 'BayesSubOpt_OnlyLambda_flat';
mdl3_projectname = 'FitAll_WT_BayesSubOpt_24July23_Final';

C = load([datadir 'ModelFits' filesep...
    mdl3_projectname filesep 'BestFit.mat']);

BestFit_mb_lambda = C.BestFit;
lambda = BestFit_mb_lambda.final_params_mean(:,5);


%% Plot

figure
histogram(lambda, NumBins = 15)

end