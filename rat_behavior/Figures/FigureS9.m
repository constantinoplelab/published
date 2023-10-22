function FigureS9(datadir, codedir)
%FigureS9 -  Inferential model identifies mistaken inferences during mixed 
% blocks across rats
%   datadir = directory of dataset
%   codedir = director of code

%% Set up paths
s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codedir, 'IgnoreCase', ispc);

if ~onPath % only add code dir to path if it already isn't
    addpath(genpath(codedir))
end

%% Load data
a = load([datadir filesep 'ratList.mat']);
ratList = a.ratList; % list of rats

mdl_projectname = 'FitAll_WT_Bayes_24July23_Final';
A = load([datadir 'ModelFits' filesep...
    mdl_projectname filesep 'BestFit.mat']);
BestFit_mb = A.BestFit;

%% Analysis

% Pre-allocate matricies
[misInfLoWTmean, misInfHiWTmean,...
    corrInfMiWTmean, corrInfHiWTmean, corrInfLoWTmean] =...
    deal(nan(length(ratList), 5));

for rr = 1:length(ratList)
    disp(rr)
    ratname = ratList{rr};

    % Pull behavioral data
    ratTrial = BestFit_mb.(ratname).Test.ratTrial;
    blk = ratTrial.block;
    rew = convertreward(ratTrial.reward);
    wt = ratTrial.wait_time;
    iscatch = ratTrial.catch;

    % Pull inferred blk from model fit
    blkInf = BestFit_mb.(ratname).Test.BlkInf;

    % Find mixed blk catch trials where the inferred blk is low
    infLo = blk==1 & blkInf==3 & iscatch;
    misInfLoWTmean(rr,:) = arrayfun(@(r)...
        mean(wt(infLo & rew==r)), 1:5); % Average WT
    
    % Find mixed blk catch trials where the inferred blk is high
    infHi = blk==1 & blkInf==2 & iscatch;
    misInfHiWTmean(rr,:) = arrayfun(@(r)...
        mean(wt(infHi & rew==r)), 1:5); % Average WT

    % Find mixed blk catch trials where the inferred blk is also mixed
    infMi = blk==1 & blkInf==1 & iscatch;
    corrInfMiWTmean(rr,:) = arrayfun(@(r)...
        mean(wt(infMi & rew==r)), 1:5);
    corrInfHiWTmean(rr,:) = arrayfun(@(r)...
        mean(wt(blk==2 & blkInf==2 & iscatch & rew==r)), 1:5);
    corrInfLoWTmean(rr,:) = arrayfun(@(r)...
        mean(wt(blk==3 & blkInf==3 & iscatch & rew==r)), 1:5);
end

p = signrank(misInfHiWTmean(:,3), misInfLoWTmean(:,3));

%% Plot
figure; hold on
l(1) = shadedErrorBar(1:5,...
    mean(misInfLoWTmean, 'omitnan'),...
    std(misInfLoWTmean, 'omitnan')./sqrt(sum(~isnan(misInfLoWTmean))),...
    lineProps='-b');
l(2) = shadedErrorBar(1:5,...
    mean(misInfHiWTmean, 'omitnan'),...
    std(misInfHiWTmean, 'omitnan')./sqrt(sum(~isnan(misInfHiWTmean))),...
    lineProps='-r');
l(3) = shadedErrorBar(1:5,...
    mean(corrInfMiWTmean, 'omitnan'),...
    std(corrInfMiWTmean, 'omitnan')./sqrt(sum(~isnan(corrInfMiWTmean))),...
    lineProps='-k');
title('Wait time in mixed blocks')

xlim([0.5 5.5])

l(1).mainLine.DisplayName = 'Inferred Low';
l(2).mainLine.DisplayName = 'Inferred High';
l(3).mainLine.DisplayName = 'Inferred Mixed';
legend(location='best')

xlabel('Wait time ratio')
title(p)

set(gcf, Position=[1121 518 560 427])
end


