function FigureS3(datadir, codedir)
%FigureS3 -  Wait times are not affected by previous trial outcome
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

%% Process data

a = load([datadir filesep 'ratList.mat']);
ratList = a.ratList; % list of rats

% block transition dynamic parameters
twin = 20; % number of trials around each block to pull
smoothfactor = 10; % size of smoothing window

% previous trial regression parameters
nback = 10; % Number of previous trials to include

[WTbyRew_RewardedMi, WTbyRew_RewardedHi, WTbyRew_RewardedLo,...
    WTbyRew_UnrewardedMi, WTbyRew_UnrewardedHi, WTbyRew_UnrewardedLo,...
    WTbyRew_AllMi, WTbyRew_AllHi, WTbyRew_AllLo] =...
    deal(nan(length(ratList), 5));

% Pre-allocation matricies
[ltomRew, htomRew, mtolRew, mtohRew,...
    ltomUnrew, htomUnrew, mtolUnrew, mtohUnrew,...
    ltomAll, htomAll, mtolAll, mtohAll] =...
    deal(nan(length(ratList), 2*twin+1));

[betasRew, betasUnrew, betasAll] = deal(nan(length(ratList), nback+1));

for rr = 1:length(ratList)
    fprintf('%d out of %d\n', rr, length(ratList));

    % Load Data
    fname = strcat(['ratTrial_', ratList{rr}, '.mat']);
    a = load([datadir 'A_Structs_Final' filesep fname]); A = a.A;

% Pull behavioral data
    wt = A.wait_time; % wait time
    postrew = [0; A.hits(1:end-1)]; % post rewarded trials (Hit = 1)
    postunrew = [0; ~A.hits(1:end-1)]; % post unrewarded trials (Hit = 0)
    rew = convertreward(A.reward); % reward (converted to 1:5)
    iscatch = A.catch & A.optout; % pull catch trials where rat opted out

    % Pull post-unrewarded wait time 
    usetheseUrew = iscatch & postunrew;

    WTUnrew = wt(usetheseUrew);
    blkUnrew = A.block(usetheseUrew);
    rewUnrew = rew(usetheseUrew);

    % Pull post-rewarded wait time 
    usetheseRew = iscatch & postrew;

    WTRew = wt(usetheseRew);
    blkRew = A.block(usetheseRew);
    rewRew = rew(usetheseRew);

    % Pull all wait times
    usetheseAll = iscatch;

    WTAll = wt(usetheseAll);
    blkAll = A.block(usetheseAll);
    rewAll = rew(usetheseAll);

% wait time by block
    % post-unrewarded wait time 
    WTbyRew_UnrewardedMi(rr,:) = arrayfun(@(r)...
        mean(WTUnrew(rewUnrew == r & blkUnrew == 1), 'omitnan'), 1:5);
    WTbyRew_UnrewardedHi(rr,:) = arrayfun(@(r)...
        mean(WTUnrew(rewUnrew == r & blkUnrew == 2), 'omitnan'), 1:5);
    WTbyRew_UnrewardedLo(rr,:) = arrayfun(@(r)...
        mean(WTUnrew(rewUnrew == r & blkUnrew == 3), 'omitnan'), 1:5);

    % post-rewarded wait time 
    WTbyRew_RewardedMi(rr,:) = arrayfun(@(r)...
        mean(WTRew(rewRew == r & blkRew == 1), 'omitnan'), 1:5);
    WTbyRew_RewardedHi(rr,:) = arrayfun(@(r)...
        mean(WTRew(rewRew == r & blkRew == 2), 'omitnan'), 1:5);
    WTbyRew_RewardedLo(rr,:) = arrayfun(@(r)...
        mean(WTRew(rewRew == r & blkRew == 3), 'omitnan'), 1:5);

    % all wait times
    WTbyRew_AllMi(rr,:) = arrayfun(@(r)...
        mean(WTAll(rewAll == r & blkAll == 1), 'omitnan'), 1:5);
    WTbyRew_AllHi(rr,:) = arrayfun(@(r)...
        mean(WTAll(rewAll == r & blkAll == 2), 'omitnan'), 1:5);
    WTbyRew_AllLo(rr,:) = arrayfun(@(r)...
        mean(WTAll(rewAll == r & blkAll == 3), 'omitnan'), 1:5);

% block transition dynamics
    % post-unrewarded wait time 
    [ltomUnrew(rr,:), htomUnrew(rr,:),...
        mtolUnrew(rr,:), mtohUnrew(rr,:)] =...
        block_dynamics_wt(A, twin, smoothfactor, usetheseUrew);

    % post-rewarded wait time 
    [ltomRew(rr,:), htomRew(rr,:), mtolRew(rr,:), mtohRew(rr,:)] =...
        block_dynamics_wt(A, twin, smoothfactor, usetheseRew);

    % all wait times
    [ltomAll(rr,:), htomAll(rr,:),...
        mtolAll(rr,:), mtohAll(rr,:)] =...
        block_dynamics_wt(A, twin, smoothfactor);

% previous reward regression
    % post-unrewarded wait time 
    betasUnrew(rr,:) = regress_wt_vs_rew(A, nback, false, usetheseUrew);

    % post-rewarded wait time 
    betasRew(rr,:) = regress_wt_vs_rew(A, nback, false, usetheseRew);

    % all wait times
    betasAll(rr,:) = regress_wt_vs_rew(A, nback, false);
end

% wait ratio
WTRatioUnrew = WTbyRew_UnrewardedHi(:,3)./WTbyRew_UnrewardedLo(:,3);
WTRatioRew = WTbyRew_RewardedHi(:,3)./WTbyRew_RewardedLo(:,3);
WTRatioAll = WTbyRew_AllHi(:,3)./WTbyRew_AllLo(:,3);

%% Plot 
ht = 4;
wd = 3;

figure

%--------------------------------------------------------------------------
% Wait time curves
%--------------------------------------------------------------------------
subplot(ht, wd, 2); hold on
shadedErrorBar(1:5, mean(WTbyRew_UnrewardedMi),...
    std(WTbyRew_UnrewardedMi)./sqrt(length(ratList)), lineprops={'k'})
shadedErrorBar(1:5, mean(WTbyRew_UnrewardedHi),...
    std(WTbyRew_UnrewardedMi)./sqrt(length(ratList)), lineprops={'r'})
shadedErrorBar(1:5, mean(WTbyRew_UnrewardedLo),...
    std(WTbyRew_UnrewardedMi)./sqrt(length(ratList)), lineprops={'b'})

xlim([0.5 5.5])
xticks(1:5)
xticklabels({5, 10, 20, 40, 80})
ylim([10 14.5])

title('unrewarded')

subplot(ht, wd, 1); hold on
shadedErrorBar(1:5, mean(WTbyRew_RewardedMi),...
    std(WTbyRew_RewardedMi)./sqrt(length(ratList)), lineprops={'k'})
shadedErrorBar(1:5, mean(WTbyRew_RewardedHi),...
    std(WTbyRew_RewardedMi)./sqrt(length(ratList)), lineprops={'r'})
shadedErrorBar(1:5, mean(WTbyRew_RewardedLo),...
    std(WTbyRew_RewardedMi)./sqrt(length(ratList)), lineprops={'b'})

xlim([0.5 5.5])
xticks(1:5)
xticklabels({5, 10, 20, 40, 80})
ylim([10 14.5])

title('rewarded')

subplot(ht, wd, 3); hold on
shadedErrorBar(1:5, mean(WTbyRew_AllMi),...
    std(WTbyRew_AllMi)./sqrt(length(ratList)), lineprops={'k'})
shadedErrorBar(1:5, mean(WTbyRew_AllHi),...
    std(WTbyRew_AllHi)./sqrt(length(ratList)), lineprops={'r'})
shadedErrorBar(1:5, mean(WTbyRew_AllLo),...
    std(WTbyRew_AllLo)./sqrt(length(ratList)), lineprops={'b'})

xlim([0.5 5.5])
xticks(1:5)
xticklabels({5, 10, 20, 40, 80})
ylim([10 14.5])

title('All')

%--------------------------------------------------------------------------
% wt ratio
%--------------------------------------------------------------------------
subplot(ht, wd, 5);
histogram(WTRatioUnrew)

xline(1, 'k--')
xline(mean(WTRatioUnrew), 'r')

subplot(ht, wd, 4);
histogram(WTRatioRew)

xline(1, 'k--')
xline(mean(WTRatioRew), 'r')

subplot(ht, wd, 6);
histogram(WTRatioAll)

xline(1, 'k--')
xline(mean(WTRatioAll), 'r')

%--------------------------------------------------------------------------
% Dynamics
%--------------------------------------------------------------------------
yl = [-0.2 0.3];

subplot(ht, wd, 8); hold on

fill([0 twin twin 0], [yl(1) yl(1) yl(2) yl(2)], 'k',...
    facealpha=0.1, linestyle='none')
shadedErrorBar(-twin:twin, mean(htomUnrew),...
    std(htomUnrew)./sqrt(length(ratList)),...
    lineprops={'r'})
shadedErrorBar(-twin:twin, mean(ltomUnrew),...
    std(ltomUnrew)./sqrt(length(ratList)),...
    lineprops={'b'})

ylim(yl)

subplot(ht, wd, 7); hold on
fill([0 twin twin 0], [yl(1) yl(1) yl(2) yl(2)], 'k',...
    facealpha=0.1, linestyle='none')

shadedErrorBar(-twin:twin, mean(htomRew),...
    std(htomRew)./sqrt(length(ratList)),...
    lineprops={'r'})
shadedErrorBar(-twin:twin, mean(ltomRew),...
    std(ltomRew)./sqrt(length(ratList)),...
    lineprops={'b'})

ylim(yl)

subplot(ht, wd, 9); hold on
fill([0 twin twin 0], [yl(1) yl(1) yl(2) yl(2)], 'k',...
    facealpha=0.1, linestyle='none')

shadedErrorBar(-twin:twin, mean(htomAll),...
    std(htomAll)./sqrt(length(ratList)),...
    lineprops={'r'})
shadedErrorBar(-twin:twin, mean(ltomAll),...
    std(ltomAll)./sqrt(length(ratList)),...
    lineprops={'b'})

ylim(yl)

%--------------------------------------------------------------------------
% Betas
%--------------------------------------------------------------------------
yl = [-0.01 0.22];

subplot(ht, wd, 11)
shadedErrorBar(0:nback-1,...
    mean(betasRew(:, 1:end-1)),...
    std(betasRew(:, 1:end-1))./sqrt(length(ratList)))

xlim([-0.5 nback-0.5])
ylim(yl)

subplot(ht, wd, 10)
shadedErrorBar(0:nback-1,...
    mean(betasUnrew(:, 1:end-1)),...
    std(betasUnrew(:, 1:end-1))./sqrt(length(ratList)))

xlim([-0.5 nback-0.5])
ylim(yl)

subplot(ht, wd, 12)
shadedErrorBar(0:nback-1,...
    mean(betasAll(:, 1:end-1)),...
    std(betasAll(:, 1:end-1))./sqrt(length(ratList)))

xlim([-0.5 nback-0.5])
ylim(yl)

set(gcf, Position=[749 80 932 865],...
    renderer='painters')

shg
end