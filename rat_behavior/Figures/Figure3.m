function [BestFit_mf, BestFit_mb,...
    BestFit_mb_lambda] = Figure3(datadir, codedir)
%Figure3 - Computational modeling reveals distinct value computations for 
% wait time and trial initiation
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

%% Load model fits - takes a while
disp('loading data. may take a while.')
mdl1 = 'AlphaPH';
mdl1_projectname = 'FitAll_WT_AlphaPH_24July23_Final';

mdl2 = 'Bayes';
mdl2_projectname = 'FitAll_WT_Bayes_24July23_Final';

mdl3 = 'BayesSubOpt_OnlyLambda_flat';
mdl3_projectname = 'FitAll_WT_BayesSubOpt_24July23_Final';

A = load([datadir 'ModelFits' filesep...
    mdl1_projectname filesep 'BestFit.mat']);
B = load([datadir 'ModelFits' filesep...
    mdl2_projectname filesep 'BestFit.mat']);
C = load([datadir 'ModelFits' filesep...
    mdl3_projectname filesep 'BestFit.mat']);

BestFit_mf = A.BestFit;
BestFit_mb = B.BestFit;
BestFit_mb_lambda = C.BestFit;

a = load([datadir filesep 'ratList.mat']);
ratList = a.ratList; % list of rats

%% Generate example kappas from each model
% generate example kappas with handchosen parameters
[Kappa_MF, Kappa_MB] = generateExampleKappas(A);

%% Wait time by catch probability

probCatchVec = [0.1, 0.15, 0.2, 0.25, 0.35]; % Possible catch probabilities

% Model simulation parameters
tau = 2.5; % exponential distribution mean
D = 1; % Scale parameter (hand-chosen)
kappa = 0.1; % value of the environment (hand-chosen)

% Preallocate matricies
WTByProbCatch = cell(length(probCatchVec), 1);
[WTByProbCatch{:}] = deal(nan(length(ratList), 5));

wtByPCatchMdl = nan(length(probCatchVec), 5);

% Loop over rats
for rr = 1:length(ratList)
    fprintf('%d out of %d\n', rr, length(ratList));

    % Load Data
    a = load([datadir 'A_Structs_Final' filesep...
        'ratTrial_' ratList{rr} '.mat']);
    A = a.A;

    % Pull reward
    rew = convertreward(A.reward); % Convert to 1:5

    % Drop first 1000 trials- first trials in stage 8 have higher catch
    % probability but rats are still learning, so drop these trials
    dontdrop = true(size(A.reward));
    dontdrop(1:1000) = false;

    % Loop over catch probabilities
    for pp = 1:length(probCatchVec)
        % Find mixed block catch trials with the catch probability
        usethese = A.block == 1 &...
            A.prob_catch==probCatchVec(pp) &...
            A.catch &...
            dontdrop;

        % Average wait times as a function of reward
        WTByProbCatch{pp}(rr,:) =...
            arrayfun(@(r)...
            mean(A.wait_time(usethese & rew==r), 'omitnan'),...
            1:5);
    end
end

% Average across rats for each catch probability
wtByPCatchRat =...
    cell2mat(cellfun(@(p) mean(p, 'omitnan'), WTByProbCatch,...
    'UniformOutput', false));

% Simulate model
for pp = 1:length(probCatchVec)
    pcatch = probCatchVec(pp);
    C = 1-pcatch;

    % Wait time from Equation (1)
    wtByPCatchMdl(pp,:) = arrayfun(@(r)...
        D*tau*log(C/(1-C)*(r-kappa*tau)/(kappa*tau)),...
        [5 10 20 40 80]);
end

%% Model Comparsion

% model comparison (BIC = Bayes Information Criteria)
[~, ~, deltaBIC, pvals] =...
    model_comparison(mdl2_projectname, mdl2,...
    mdl1_projectname, mdl1,...
    true, true, true, false, BestFit_mb, BestFit_mf);

%% Generating Model Transitions

% block transition dynamic parameters - wt model
twin1 = 20; % number of trials around each block to pull
smoothfactor1 = 10; % size of smoothing window

% block transition dynamic parameters - latency model
twin2 = 20; % number of trials around each block to pull
smoothfactor2 = 15; % size of smoothing window

% Number of model simulations
Nsims = 50;

% Parameters for the simulations
params_mf = [0.1042, 0.2527]; % retrospective modle
params_mb = [0.5, 0.9, 0.3 0.9]; % inferential model

% Use the behavioral data from J004 for simluations
[A_mf, A_mb] = deal(BestFit_mb.J004.All.ratTrial);

% Preallocate matricies
[ltom_mf, htom_mf] = deal(nan(Nsims, 2*twin2+1));
[ltom_mb, htom_mb] = deal(nan(Nsims, 2*twin1+1));

rng(724) % Set random seed for reproducibility 

% Simulate
for rr = 1:Nsims
    fprintf('%d out of %d\n', rr, length(Nsims));

    % Simulate inferential model
    [~, A_mb.wait_time] = GenerateSynthData_Bayes(params_mb, A_mb,...
        'logn', true, 8);
    [ltom_mb(rr,:), htom_mb(rr,:)] =...
        block_dynamics_wt(A_mb, twin1, smoothfactor1);
    
    % Simulate retrospective model
    [~, A_mf.ITI] = GenerateLatencyData_AlphaPHModel(params_mf, A_mf,...
        true, 'logn', 8);
    [ltom_mf(rr,:), htom_mf(rr,:)] =...
        block_dynamics_latency(A_mf, twin2, smoothfactor2);
end

%% Model vs. rat data - wt

% Example rat for just mixed blocks
ratname1 = 'J069'; % hand-chosen based on performance
[blk_rat_mean1, blk_rat_sem1, blk_mdl_mean1, blk_mdl_sem1] =...
    ratAndMdlWTCurves(ratname1, BestFit_mb);

% Example rat for just high and low blocks
ratname2 = 'G051'; % hand-chosen based on performance
[blk_rat_mean2, blk_rat_sem2, blk_mdl_mean2, blk_mdl_sem2] =...
    ratAndMdlWTCurves(ratname2, BestFit_mb);

%% Model ti time

% Parameters for model simulation (hand-chosen)
params_mf_ex = [0.3, 0.62]; 

% Use J031 data for model simulation
A = BestFit_mb.J031.All.ratTrial;

% Pull block
blk = A.block;
usethese_hi = blk==2; % find high blocks
usethese_lo = blk==3; % find low blocks

% generate model trial initation times
[~, L_mdl] =...
    GenerateLatencyData_AlphaPHModel(params_mf_ex, A, true, 'logn', 1.5);

% Average and SEM over blocks
lat_mdl =...
    [mean(L_mdl(usethese_lo), 'omitnan'),...
    mean(L_mdl(usethese_hi), 'omitnan');...
    std(L_mdl(usethese_lo), 'omitnan')./sqrt(sum(usethese_lo)),...
    std(L_mdl(usethese_hi), 'omitnan')./sqrt(sum(usethese_hi))];

%% Model conditional wait time

% Preallocate matricies
[cond_wait_time_mean, cond_iti_mean] = deal(nan(length(ratList), 2));

% Loop over rats
for rr = 1:length(ratList)
    fprintf('%d out of %d\n', rr, length(ratList));
    ratname = ratList{rr};

    % Pull ratTrial for rat
    A = BestFit_mb.(ratname).All.ratTrial;

    % Set wait time equal to model predicition wait time (function for
    % conditional wait/trial initiation time use A struct)
    A.wait_time = BestFit_mb.(ratname).All.WTMdl;

    % Set wait time equal to model predicition wait time (function for
    % conditional wait/trial initiation time use A struct)
    [~, A.ITI] =...
        GenerateLatencyData_AlphaPHModel(...
        BestFit_mf.(ratname).final_params(1,:), A, true,...
        'logn', 1.5);

    % Calc conditional wait/trial initiaiton times (same function as rats)
    cond_wait_time_mean(rr,:) = waittime_20ul_by_previous_vol(A);
    cond_iti_mean(rr,:) = latency_post20ul_by_previous_vol(A);
end

%% Block dynamics by lambda

% Shuffle test parameters
N = 60; % Percentile to split high/low lambda
Nshuffle = 500; % Number of shuffles to perform

% block transition dynamic parameters
twin3 = 30; % number of trials around each block to pull
smoothfactor = 10; % size of smoothing window

% Preallocate matricies
[ltomWT, htomWT, mtolWT, mtohWT,...
    ltomTI, htomTI, mtolTI, mtohTI] =...
    deal(nan(length(ratList), 2*twin3+1));
lambda = nan(size(ratList));

% Loop over all rats
for rr = 1:length(ratList)
    fprintf('%d out of %d\n', rr, length(ratList));

    % Pull rat data
    A = BestFit_mb_lambda.(ratList{rr}).All.ratTrial;

    % Pull block dynamics for wait time and trial initiation time
    [ltomWT(rr,:), htomWT(rr,:), mtolWT(rr,:), mtohWT(rr,:)] =...
        block_dynamics_wt(A, twin3, smoothfactor);
    [ltomTI(rr,:), htomTI(rr,:), mtolTI(rr,:), mtohTI(rr,:)] =...
        block_dynamics_latency(A, twin3, smoothfactor);

    % Pull that rat's lambda
    lambda(rr) = BestFit_mb_lambda.(ratList{rr}).final_params(1,5);
end

% Find high/low lambda rats
hiLambda = lambda > prctile(lambda, N); % high lambdas (good inference)
loLambda = lambda < prctile(lambda, 100-N); % low lambdas (poor inference)

% Pull data for high and low lambda rats
mtolWTlo = mtolWT(loLambda,:); mtolWThi = mtolWT(hiLambda,:);
mtohWTlo = mtohWT(loLambda,:); mtohWThi = mtohWT(hiLambda,:);

mtolTIlo = mtolTI(loLambda,:); mtohTIhi = mtohTI(loLambda,:);
mtolTIhi = mtolTI(hiLambda,:); mtohTIlo = mtohTI(hiLambda,:);

% Perform shuffle test
[~, ~, pvalsWTmtol] =...
    myLambdaShuffleTest(mtolWTlo, mtolWThi, Nshuffle, twin3);
[~, ~, pvalsWTmtoh] =...
    myLambdaShuffleTest(mtohWTlo, mtohWThi, Nshuffle, twin3);
[~, ~, pvalsTImtol] =...
    myLambdaShuffleTest(mtolTIlo, mtolTIhi, Nshuffle, twin3);
[~, ~, pvalsTImtoh] =...
    myLambdaShuffleTest(mtohTIlo, mtohTIhi, Nshuffle, twin3);

% Print p-values 
disp([pvalsWTmtol; pvalsWTmtoh; pvalsTImtol; pvalsTImtoh])

%% Plot
clear l

mbColor = '#a16ae8';
mfColor = '#50C878';
cmap = generate_color_gradient([0.85 0.85 0.85], [0 0 0],...
    length(probCatchVec));

ht = 4;
wd = 4;

norm = @(x) x./max(abs(x));
myNorm = @(x) x-mean(x(1:twin3));

figure; colororder(cmap)
%--------------------------------------------------------------------------
% Model eq
%--------------------------------------------------------------------------
subplot(ht, wd, 1); hold on
title('Model eq')

%--------------------------------------------------------------------------
% Wait time in mixed blocks
%--------------------------------------------------------------------------
subplot(ht, wd, 2); hold on
shadedErrorBar(1:5, blk_rat_mean1(1,:), blk_rat_sem1(1,:),...
    lineprops={'color', '#980036'})
shadedErrorBar(1:5, blk_mdl_mean1(1,:), blk_mdl_sem1(1,:),...
    lineprops={'color', '#ff8b00'})

title('Model in Mixed Blocks')
xlim([0.5 5.5])
%--------------------------------------------------------------------------
% Wait time by pcatch - rat
%--------------------------------------------------------------------------
subplot(ht, wd, 3); hold on
plot(1:5, wtByPCatchRat')

title('rat by pcatch')
xlim([0.5 5.5])
ylim([8.75 18])
%--------------------------------------------------------------------------
% Wait time by pcatch - mdl
%--------------------------------------------------------------------------
subplot(ht, wd, 4); hold on
plot(1:5, wtByPCatchMdl')
xlim([0.5 5.5])
ylim([8 20.5])

%--------------------------------------------------------------------------
% Bayes kappa
%--------------------------------------------------------------------------

subplot(ht, wd, 5); hold on
yl = [0.23 0.9];

fill([0 40 40 0], [yl(1) yl(1) yl(2) yl(2)],...
    'k', facealpha=0.15, edgecolor='none')
fill([40 80 80 40], [yl(1) yl(1) yl(2) yl(2)],...
    'b', facealpha=0.15, edgecolor='none')
fill([80 120 120 80], [yl(1) yl(1) yl(2) yl(2)],...
    'k', facealpha=0.15, edgecolor='none')
fill([120 160 160 120], [yl(1) yl(1) yl(2) yl(2)],...
    'r', facealpha=0.15, edgecolor='none')
fill([160 200 200 160], [yl(1) yl(1) yl(2) yl(2)],...
    'k', facealpha=0.15, edgecolor='none')

plot(Kappa_MB, '-', color=mbColor, linewidth=2.5,...
    markerfacecolor=mbColor, markeredgecolor='none', markersize=4)

ylim(yl)

yticks([])
xticks([])

xlabel('Trial')
ylabel('Opportunity Cost')

%--------------------------------------------------------------------------
% WT Model Dynamics
%--------------------------------------------------------------------------
twin = 20;
yl = [-1.22 1.22];

subplot(ht, wd, 6); hold on
fill([0 twin twin 0], [yl(1) yl(1) yl(2) yl(2)],...
    'k', facealpha=0.1, linestyle='none')

l(1) = shadedErrorBar(-twin1:twin1, norm(mean(ltom_mb)),...
    std(ltom_mb), lineprops={'b'});
l(2) = shadedErrorBar(-twin1:twin1, norm(mean(htom_mb)),...
    std(ltom_mb), lineprops={'r'});

xlim([-20 20])
ylim(yl)

ylabel('Wait-time (z-score)')

%--------------------------------------------------------------------------
% Model vs. Rat High/Low Comparison
%--------------------------------------------------------------------------

subplot(ht, wd, 7)
l(5) = shadedErrorBar(1:5, blk_rat_mean2(2,:), blk_rat_sem2(2,:),...
    lineprops={'color', '#980036'});
l(6) = shadedErrorBar(1:5, blk_rat_mean2(3,:), blk_rat_sem2(3,:),...
    lineprops={'color', '#980036'});

xlim([0.5 5.5])
l(7) = shadedErrorBar(1:5, blk_mdl_mean2(2,:), blk_mdl_sem2(2,:),...
    lineprops={'color', '#ff8b00'});
l(8) = shadedErrorBar(1:5, blk_mdl_mean2(3,:), blk_mdl_sem2(3,:),...
    lineprops={'color', '#ff8b00'});

set(l(5).mainLine, DisplayName='rat')
set(l(7).mainLine, DisplayName='model')
legend([l(5).mainLine, l(7).mainLine], location='best')
title(ratname2)

xticks(1:5)

ylabel('Wait time (s)')
title(ratname)

%--------------------------------------------------------------------------
% Conditional wt model
%--------------------------------------------------------------------------
subplot(ht, wd, 8); hold on
plot(1:2, cond_wait_time_mean, color=[0.7 0.7 0.7])
plot(1:2, mean(cond_wait_time_mean), 'k')

xticks(1:2)
xticklabels({'<20', '>20'})
xlim([0.5 2.5])
ylim([-1 1])
title('model conditional wt')
%--------------------------------------------------------------------------
% MF kappa
%--------------------------------------------------------------------------
subplot(ht, wd, 9); hold on
yl = [0.34 1.05];

fill([0 40 40 0], [yl(1) yl(1) yl(2) yl(2)],...
    'k', facealpha=0.15, edgecolor='none')
fill([40 80 80 40], [yl(1) yl(1) yl(2) yl(2)],...
    'b', facealpha=0.15, edgecolor='none')
fill([80 120 120 80], [yl(1) yl(1) yl(2) yl(2)],...
    'k', facealpha=0.15, edgecolor='none')
fill([120 160 160 120], [yl(1) yl(1) yl(2) yl(2)],...
    'r', facealpha=0.15, edgecolor='none')
fill([160 200 200 160], [yl(1) yl(1) yl(2) yl(2)],...
    'k', facealpha=0.15, edgecolor='none')

plot(Kappa_MF, '-', color=mfColor, linewidth=2.5,...
    markerfacecolor=mfColor, markeredgecolor='none', markersize=4)

plot([81 120], [mean(Kappa_MF(81:120)) mean(Kappa_MF(81:120))],...
    'k--')
plot([161 200], [mean(Kappa_MF(161:200)) mean(Kappa_MF(161:200))],...
    'k--')

% ylim(yl)

yticks([])
xticks([])

xlabel('Trial')
ylabel('Opportunity Cost')

%--------------------------------------------------------------------------
% Lat Model Dynamics
%--------------------------------------------------------------------------

subplot(ht, wd, 10); hold on
yl = [-1.075 1.075];
fill([0 twin twin 0], [yl(1) yl(1) yl(2) yl(2)],...
    'k', facealpha=0.1, linestyle='none')

l(3) = shadedErrorBar(-twin2:twin2, norm(mean(ltom_mf)),...
    std(ltom_mf)./sqrt(size(ltom_mf, 1)), lineprops={'b'});
l(4) = shadedErrorBar(-twin2:twin2, norm(mean(htom_mf)),...
    std(htom_mf)./sqrt(size(htom_mf, 1)), lineprops={'r'});

xlim([-20 20])
ylim(yl)

ylabel('Trial initiation time (z-score)')

%--------------------------------------------------------------------------
% Model vs. Rat High/Low Comparison
%--------------------------------------------------------------------------
subplot(ht, wd, 11); hold on
% errorbar([1 2]-0.015, lat_rat(1,:), lat_rat(2,:),...
%     'ko', markerfacecolor='k', linewidth=2, capsize=0)
errorbar(1, lat_mdl(1,1), lat_mdl(2,1),...
    'bo', markerfacecolor='b', linewidth=2, capsize=0)
errorbar(2, lat_mdl(1,2), lat_mdl(2,2),...
    'ro', markerfacecolor='r', linewidth=2, capsize=0)

ylim([2.5 5])
xlim([0.5 2.5])

%--------------------------------------------------------------------------
% By previous volume, model?
%--------------------------------------------------------------------------
subplot(ht, wd, 12); hold on
plot(1:2, cond_iti_mean, color=[0.7 0.7 0.7])
plot(1:2, mean(cond_iti_mean), 'k')

xticks(1:2)
xticklabels({'<20', '>20'})
xlim([0.5 2.5])
ylim([-0.11 0.11])
title('model conditional iti')

%--------------------------------------------------------------------------
% Model Comparison
%--------------------------------------------------------------------------

nbins = 22;
xs1 = linspace(0, 1000, nbins/2);
xs2 = linspace(-1000, 0, nbins/2);
binedges = [xs2(1:end-1), xs1];

subplot(ht, wd, 13)
h = histogram(deltaBIC, BinEdges=binedges);

xline(0, 'k--', linewidth=2, alpha=1)
xline(mean(deltaBIC, 'omitnan'), 'r--', linewidth=2, alpha=1)
title(['\Delta BIC ' num2str(pvals(3))])

set(h, BinLimits = [-max(abs(h.BinLimits)), max(abs(h.BinLimits))])

xlabel('mf-mb')
ylabel('n (rats)')

%--------------------------------------------------------------------------
% Lambda percentile dynamics, wt
%--------------------------------------------------------------------------

yl_wt = [-0.1 0.19];

subplot(ht, wd, 15); hold on
text(-10, 0.15, {pvalsWTmtol, pvalsWTmtoh})

fill([-twin 0 0 -twin], [yl_wt(1) yl_wt(1) yl_wt(2) yl_wt(2)], 'k',...
    facealpha=0.1, linestyle='none')

l(9) = shadedErrorBar(-twin3:twin3, myNorm(mean(mtolWTlo)),...
    std(mtolWTlo)./sqrt(size(mtolWTlo, 1)),...
    lineprops={'b'});
l(10) = shadedErrorBar(-twin3:twin3, myNorm(mean(mtolWThi)),...
    std(mtolWThi)./sqrt(size(mtolWThi, 1)),...
    lineprops={'b--'});
l(11) = shadedErrorBar(-twin3:twin3, myNorm(mean(mtohWTlo)),...
    std(mtohWTlo)./sqrt(size(mtohWTlo, 1)),...
    lineprops={'r'});
l(12) = shadedErrorBar(-twin3:twin3, myNorm(mean(mtohWThi)),...
    std(mtohWThi)./sqrt(size(mtohWThi, 1)),...
    lineprops={'r--'});

xlim([-15 25])
ylim(yl_wt)


%--------------------------------------------------------------------------
% Lambda percentile dynamics, ti time
%--------------------------------------------------------------------------
yl_ti = [-0.075 0.2];

subplot(ht, wd, 16); hold on
text(-10, 0.15, {pvalsTImtol, pvalsTImtoh})

fill([-twin 0 0 -twin], [yl_ti(1) yl_ti(1) yl_ti(2) yl_ti(2)], 'k',...
    facealpha=0.1, linestyle='none')

l(13) = shadedErrorBar(-twin3:twin3, myNorm(mean(mtolTIlo)),...
    std(mtolTIlo)./sqrt(size(mtolTIlo, 1)),...
    lineprops={'b'});
l(14) = shadedErrorBar(-twin3:twin3, myNorm(mean(mtolTIhi)),...
    std(mtolTIhi)./sqrt(size(mtolTIhi, 1)),...
    lineprops={'b--'});
l(15) = shadedErrorBar(-twin3:twin3, myNorm(mean(mtohTIlo)),...
    std(mtohTIlo)./sqrt(size(mtohTIlo, 1)),...
    lineprops={'r'});
l(16) = shadedErrorBar(-twin3:twin3, myNorm(mean(mtohTIhi)),...
    std(mtohTIhi)./sqrt(size(mtohTIhi, 1)),...
    lineprops={'r--'});

xlim([-15 20])
ylim(yl_ti)

set(gcf, Position=[664 117 1017 830], renderer='painters');
arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), l)



end