function FigureS11(datadir, codedir)
%FigureS11 -  Differential wait time dynamics based on lambda from 
% sub-optimal Bayes model are robust across a range of percentiles.
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
a = load([datadir filesep 'ratList.mat']);
ratList = a.ratList; % list of rats

mdl3_projectname = 'FitAll_WT_BayesSubOpt_24July23_Final';
C = load([datadir 'ModelFits' filesep...
    mdl3_projectname filesep 'BestFit.mat']);
BestFit_mb_lambda = C.BestFit;

%% Analyze data
% Shuffle test parameters
Nshuffle = 500; % number of shuffles

% block transition dynamic parameters
twin = 30; % number of trials around each block to pull
smoothfactor = 10; % size of smoothing window

% Pre-allocate matricies
[ltomWT, htomWT, mtolWT, mtohWT] = deal(nan(length(ratList), 2*twin+1));
lambda = nan(size(ratList));

% Loop over rats
for rr = 1:length(ratList)
    fprintf('%d out of %d\n', rr, length(ratList));
    
    A = BestFit_mb_lambda.(ratList{rr}).All.ratTrial;

    [ltomWT(rr,:), htomWT(rr,:), mtolWT(rr,:), mtohWT(rr,:)] =...
        block_dynamics_wt(A, twin, smoothfactor);

    lambda(rr) = BestFit_mb_lambda.(ratList{rr}).final_params(1,5);
end

%%
myNorm = @(x) x-mean(x(1:twin));

Ns = [60, 65, 70];

[mtolWTlo, mtolWThi, mtohWTlo, mtohWThi,...
    pvalsWTmtol, pvalsWTmtoh] =...
    deal(cell(size(Ns)));

for nn = 1:length(Ns)
    N = Ns(nn);

    hiLambda = lambda > prctile(lambda, N);
    loLambda = lambda < prctile(lambda, 100-N);

    mtolWTlo{nn} = mtolWT(loLambda,:); mtolWThi{nn} = mtolWT(hiLambda,:);
    mtohWTlo{nn} = mtohWT(loLambda,:); mtohWThi{nn} = mtohWT(hiLambda,:);

    % Shuffle test
    [~, ~, pvalsWTmtol{nn}] =...
        myLambdaShuffleTest(mtolWTlo{nn}, mtolWThi{nn}, Nshuffle, twin);

    [~, ~, pvalsWTmtoh{nn}] =...
        myLambdaShuffleTest(mtohWTlo{nn}, mtohWThi{nn}, Nshuffle, twin);
end

%% Plot

figure
for nn = 1:length(Ns)

    yl_wt = [-0.2 0.2];

    subplot(1, 3, nn); hold on
    fill([-twin 0 0 -twin], [yl_wt(1) yl_wt(1) yl_wt(2) yl_wt(2)], 'k',...
        facealpha=0.1, linestyle='none')

    shadedErrorBar(-twin:twin, myNorm(mean(mtolWTlo{nn})),...
        std(mtolWTlo{nn})./sqrt(size(mtolWTlo{nn}, 1)),...
        lineprops={'b'});
    shadedErrorBar(-twin:twin, myNorm(mean(mtolWThi{nn})),...
        std(mtolWThi{nn})./sqrt(size(mtolWThi{nn}, 1)),...
        lineprops={'b--'});
    shadedErrorBar(-twin:twin, myNorm(mean(mtohWTlo{nn})),...
        std(mtohWTlo{nn})./sqrt(size(mtohWTlo{nn}, 1)),...
        lineprops={'r'});
    shadedErrorBar(-twin:twin, myNorm(mean(mtohWThi{nn})),...
        std(mtohWThi{nn})./sqrt(size(mtohWThi{nn}, 1)),...
        lineprops={'r--'});
    
    title([num2str(Ns(nn)) ' / ' num2str(100-Ns(nn))])
    xlim([-15 25])
    ylim(yl_wt)

    disp([num2str(Ns(nn)) ' / ' num2str(100-Ns(nn))])
    disp(['mixed to low ps: ' num2str(pvalsWTmtol{nn})])
    disp(['mixed to high ps: ' num2str(pvalsWTmtoh{nn})])
end

set(gcf, Position=[398 549 1473 420])