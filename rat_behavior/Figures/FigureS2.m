function FigureS2(datadir, codedir)
%FigureS2 -  Trial initiation times depend on previous trial outcome
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

%% Supp fig - ITI for rewarded vs. unrewarded trials

a = load([datadir filesep 'ratList.mat']);
ratList = a.ratList; % list of rats

% block transition dynamic parameters
twin = 20; % number of trials around each block to pull
smoothfactor = 10; % size of smoothing window

% previous trial regression parameters
nback = 10; % Number of previous trials to include

% Pre-allocation matricies
[LatRew, LatUnrew, LatAll] = deal(nan(length(ratList), 2));
[LatRewRatio, LatUnrewRatio, LatAllRatio] = deal(nan(size(ratList)));

[ltomRew, htomRew, ltomUnrew, htomUnrew, ltomAll, htomAll] =...
    deal(nan(length(ratList), 2*twin+1));

[betasRew, betasUnrew, betasAll] = deal(nan(length(ratList), nback+1));

[rewardedITIbyPrevVol, UnrewardedITIbyPrevVol, AllITIbyPrevVol] =...
    deal(nan(length(ratList), 5));

% Process each rat's data
for rr = 1:length(ratList)
    fprintf('%d out of %d\n', rr, length(ratList));

    % Load Data
    fname = strcat(['ratTrial_', ratList{rr}, '.mat']);
    a = load([datadir 'A_Structs_Final' filesep, fname]);
    A = a.A;

% Pull behavioral data
    l = A.ITI; % latency
    postrew = [0; A.hits(1:end-1)]; % post rewarded trials (Hit = 1)
    postunrew = [0; ~A.hits(1:end-1)]; % post unrewarded trials (Hit = 0)
    rew = convertreward(A.reward); % reward (converted to 1:5)
    prevrew = [nan; rew(1:end-1)]; % previous reward (converted to 1:5) 

    % process latency
    l(l>prctile(l,99)) = NaN; % Filter large outliers 
    lz = (l-mean(l, 'omitnan'))./std(l, 'omitnan'); % z-score

    % Pull Post-Rewarded Latency
    usetheseRew = ~isnan(l) & postrew;

    LRew = l(usetheseRew);
    LRewZ = lz(usetheseRew);
    blkRew = A.block(usetheseRew);
    prevRewRew = prevrew(usetheseRew);

    % Pull Post-Unrewarded Latency
    usetheseUnrew = ~isnan(l) & postunrew;

    LUnrew = l(usetheseUnrew);
    LUnrewZ = lz(usetheseUnrew);
    blkUnrew = A.block(usetheseUnrew);
    prevRewUnrew = prevrew(usetheseUnrew);

    % Pull All Latency
    usetheseAll = ~isnan(l);

    LAll = l(usetheseAll);
    LAllZ = lz(usetheseAll);
    blkAll = A.block(usetheseAll);
    prevRewAll = prevrew(usetheseAll);

% Trial initiation time by previous reward
    % Post-rewarded trials
    rewardedITIbyPrevVol(rr,:) =...
        arrayfun(@(r)...
        mean(LRewZ(prevRewRew == r & blkRew==1), 'omitnan'), 1:5);
    
    % Post-unrewarded trials
    UnrewardedITIbyPrevVol(rr,:) =...
        arrayfun(@(r)...
        mean(LUnrewZ(prevRewUnrew == r & blkUnrew==1), 'omitnan'), 1:5);

    % All trials
    AllITIbyPrevVol(rr,:) =...
        arrayfun(@(r)...
        mean(LAllZ(prevRewAll == r & blkAll==1), 'omitnan'), 1:5);

% Trial initiation time by block
    % Post-rewarded trials
    LatRew(rr,:) =...
        [mean(LRewZ(blkRew==3), 'omitnan'),...
        mean(LRewZ(blkRew==2), 'omitnan')];

    % Post-unrewarded trials
    LatUnrew(rr,:) =...
        [mean(LUnrewZ(blkUnrew==3), 'omitnan'),...
        mean(LUnrewZ(blkUnrew==2), 'omitnan')];
    
    % All trials
    LatAll(rr,:) =...
        [mean(LAllZ(blkAll==3), 'omitnan'),...
        mean(LAllZ(blkAll==2), 'omitnan')];

% Trial initiation time ratio
    % Post-rewarded trials
    LatRewRatio(rr) =...
        mean(LRew(blkRew==2), 'omitnan')./...
        mean(LRew(blkRew==3), 'omitnan');    

    % Post-unrewarded trials
    LatUnrewRatio(rr) =...
        mean(LUnrew(blkUnrew==2), 'omitnan')./...
        mean(LUnrew(blkUnrew==3), 'omitnan');

    % All trials
    LatAllRatio(rr) =...
        mean(LAll(blkAll==2), 'omitnan')./...
        mean(LAll(blkAll==3), 'omitnan');

% Block transition dynamics
    % Post-rewarded trials
    [ltomRew(rr,:), htomRew(rr,:)] =...
        block_dynamics_latency(A, twin, smoothfactor, usetheseRew);
    
    % Post-unrewarded trials
    [ltomUnrew(rr,:), htomUnrew(rr,:)] =...
        block_dynamics_latency(A, twin, smoothfactor, usetheseUnrew);
    
    % All trials
    [ltomAll(rr,:), htomAll(rr,:)] =...
        block_dynamics_latency(A, twin, smoothfactor);

% Previous reward regression
    % Post-rewarded trials
    betasRew(rr,:) =...
        regress_latency_vs_rew(A, nback, false, true, usetheseRew);
    
    % Post-unrewarded trials
    betasUnrew(rr,:) =...
        regress_latency_vs_rew(A, nback, false, usetheseUnrew);
    
    % All trials
    betasAll(rr,:) =...
        regress_latency_vs_rew(A, nback, false, true);
end

prew = signrank(LatRew(:,1), LatRew(:,2), tail='right');
punrew = signrank(LatUnrew(:,1), LatUnrew(:,2), tail='right');
pall = signrank(LatAll(:,1), LatAll(:,2), tail='right');

%% Plot

figure
%--------------------------------------------------------------------------
% Lat by prev volume
%--------------------------------------------------------------------------
subplot(5, 3, 1); hold on
plot(1:5, mean(rewardedITIbyPrevVol))

ylim([-0.3 0.4])
xlim([0.5 5.5])
xticks(1:5)
xticklabels({'5', '10', '20' '40', '80'})

subplot(5, 3, 2); hold on
plot(1:5, mean(UnrewardedITIbyPrevVol))

ylim([-0.3 0.4])
xlim([0.5 5.5])
xticks(1:5)
xticklabels({'5', '10', '20' '40', '80'})

subplot(5, 3, 3); hold on
plot(1:5, mean(AllITIbyPrevVol))

ylim([-0.3 0.4])
xlim([0.5 5.5])
xticks(1:5)
xticklabels({'5', '10', '20' '40', '80'})

%--------------------------------------------------------------------------
% Lat by block
%--------------------------------------------------------------------------
subplot(5, 3, 4); hold on
plot([1 2], LatRew, color=[0 0 0 0.15], LineWidth=0.25)
plot([1 2], mean(LatRew), 'r', linewidth=2)

xlim([0.5 2.5])

xticks(1:2)
% xticklabels({'Low', 'High'})
title({'Post-rewarded', prew})

subplot(5, 3, 5); hold on
plot([1 2], LatUnrew, color=[0 0 0 0.15], LineWidth=0.25)
plot([1 2], mean(LatUnrew), 'r', linewidth=2)

xlim([0.5 2.5])

xticks(1:2)
xticklabels({'Low', 'High'})
title({'Post-unrewarded', punrew})

subplot(5, 3, 6); hold on
plot([1 2], LatAll, color=[0 0 0 0.15], LineWidth=0.25)
plot([1 2], mean(LatAll), 'r', linewidth=2)

xlim([0.5 2.5])

xticks(1:2)
xticklabels({'Low', 'High'})
title({'All', pall})

%--------------------------------------------------------------------------
% Block ratio
%--------------------------------------------------------------------------
subplot(5, 3, 7); hold on
histogram(LatRewRatio)

xline(1, 'k--', linewidth=2)
xlabel('ITI Ratio')

subplot(5, 3, 8); hold on
histogram(LatUnrewRatio)

xline(1, 'k--', linewidth=2)
xlabel('ITI Ratio')

subplot(5, 3, 9); hold on
histogram(LatAllRatio)

xline(1, 'k--', linewidth=2)
xlabel('ITI Ratio')

%--------------------------------------------------------------------------
% Block Transition
%--------------------------------------------------------------------------
subplot(5, 3, 10); hold on
plot(-twin:twin, mean(ltomRew), 'b', linewidth=2)
plot(-twin:twin, mean(htomRew), 'r', linewidth=2)

yl = [-0.2000 0.4234];

fill([0 twin twin 0], [yl(1) yl(1) yl(2) yl(2)], 'k',...
    facealpha=0.1, linestyle='none')
ylim(yl)

subplot(5, 3, 11); hold on
plot(-twin:twin, mean(ltomUnrew), 'b', linewidth=2)
plot(-twin:twin, mean(htomUnrew), 'r', linewidth=2)

fill([0 twin twin 0], [yl(1) yl(1) yl(2) yl(2)], 'k',...
    facealpha=0.1, linestyle='none')
ylim(yl)

subplot(5, 3, 12); hold on
plot(-twin:twin, mean(ltomAll), 'b', linewidth=2)
plot(-twin:twin, mean(htomAll), 'r', linewidth=2)

fill([0 twin twin 0], [yl(1) yl(1) yl(2) yl(2)], 'k',...
    facealpha=0.1, linestyle='none')
ylim(yl)

%--------------------------------------------------------------------------
% Regression betas
%--------------------------------------------------------------------------
subplot(5, 3, 13); hold on
plot(mean(betasRew(:, 2:end-1)), linewidth=2)

% ylabel('Regression Coefficient')
xlabel('Trials Back')

subplot(5, 3, 14); hold on
plot(mean(betasUnrew(:, 2:end-1)), linewidth=2)
xlabel('Trials Back')

subplot(5, 3, 15); hold on
plot(mean(betasAll(:, 2:end-1)), linewidth=2)
xlabel('Trials Back')

set(gcf, Position=[994 78 687 867],...
    renderer='painters')
end