function [wtByPCatchRat, wtByPCatchMdl] = wtByPCatch(codedir,datadir)
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here

%% Set up paths
s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codedir, 'IgnoreCase', ispc);

if ~onPath
    addpath(genpath(codedir))
end
%%

a = load([datadir 'ratList.mat']);
ratList = a.ratList;

% Catch probability vecctor
probCatchVec = [0.1, 0.15, 0.2, 0.25, 0.35];

% Parameters for model simulaiton
tau = 2.5;
D = 1;
kappa = 0.1;

% Preallocate matricies
WTByProbCatch = cell(length(probCatchVec), 1);
[WTByProbCatch{:}] = deal(nan(length(ratList), 5));

wtByPCatchMdl = nan(length(probCatchVec), 5);

for rr = 1:length(ratList)
    fprintf('%d out of %d\n', rr, length(ratList));

    % Load Data
    a = load([datadir 'A_Structs_Final' filesep...
        'ratTrial_' ratList{rr} '.mat']); A = a.A;

    % Pull reward and convert to ordinal (1:5)
    rew = convertreward(A.reward);

    % Drop the first 1000 trials of stage 8 - too early in training to use
    dontdrop = true(size(A.reward));
    dontdrop(1:1000) = false;

    % Loop over possible catch probabilities 
    for pp = 1:length(probCatchVec)
        usethese = A.block == 1 &...
            A.prob_catch==probCatchVec(pp) &...
            A.catch &...
            dontdrop;

        WTByProbCatch{pp}(rr,:) =...
            arrayfun(@(r)...
            mean(A.wait_time(usethese & rew==r), 'omitnan'),...
            1:5);
    end
end

wtByPCatchRat =...
    cell2mat(cellfun(@(p) mean(p, 'omitnan'), WTByProbCatch,...
    'UniformOutput', false));

for pp = 1:length(probCatchVec)
    pcatch = probCatchVec(pp);
    C = 1-pcatch;

    wtByPCatchMdl(pp,:) = arrayfun(@(r)...
        D*tau*log(C/(1-C)*(r-kappa*tau)/(kappa*tau)),...
        [5 10 20 40 80]);
end
end