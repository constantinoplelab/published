function [WT, Lat] = Figure2(datadir, codedir)
%Figure2 - Wait and trial initiation times use distinct estimates of the 
% value of the environment. 
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

% block transition dynamic parameters
twin = 20; % number of trials around each block to pull
smoothfactor = 10; % size of smoothing window

% Fast and slow trial intiation time tau shuffle test parameters
N = 80; % Percentile to be considered slow (or 100-N for fast)
Nshuffle = 500; % number of shuffles to perform

% previous trial regression parameters
nback = 10; % Number of previous trials to include

% Preallocate matricies
wtByPrevRew = nan(length(ratList), 2);
latByPrevRew = nan(length(ratList), 2);

pWt = nan(length(ratList), 1);
pLat = nan(length(ratList), 1);

[WT.ltom, WT.htom, Lat.ltom, Lat.htom] =...
    deal(nan(length(ratList), 2*twin+1));

tau = struct;
tau.latency = nan(length(ratList), 1);
tau.wt = nan(length(ratList), 1);

betas = struct;
betas.latency = nan(length(ratList), nback+1);
betas.wt = nan(length(ratList), nback+1);

% Set up color scheme
mbColor = '#a16ae8';
mfColor = '#50C878';

wt_color = [0.6314, 0.4157, 0.9098];
lat_color = [0.3137, 0.7843, 0.4706];

firstQColor = [249, 163, 26]/255;
lastQColor = [193, 119, 176]/255;

% Process each rat's data
for rr = 1:length(ratList)
    fprintf('%d out of %d\n', rr, length(ratList));

    % Load Data
    fname = strcat(['ratTrial_', ratList{rr}, '.mat']);
    A = load([datadir 'A_Structs_Final' filesep fname]);
    A = A.A;

    % wait time on 20 uL conditioned on the previous volume in mixed blks
    [wtByPrevRew(rr,:), pWt(rr)] = waittime_20ul_by_previous_vol(A);

    % trial initiation time conditioned on the previous volume in mixed blk
    [latByPrevRew(rr,:), ~, pLat(rr)] =...
        latency_post20ul_by_previous_vol(A);

    % Wait time dynamics around block transitions
    [WT.ltom(rr,:), WT.htom(rr,:)] =...
        block_dynamics_wt(A, twin, smoothfactor);

    % trial intiation time dynamics around block transitions
    [Lat.ltom(rr,:), Lat.htom(rr,:), Lat.mtol(rr,:), Lat.mtoh(rr,:)] =...
        block_dynamics_latency(A, twin, smoothfactor);

    % regress wait time vs. previous reward reward
    [betas.wt(rr,:), tau.wt(rr)] =...
        regress_wt_vs_rew(A, nback);

    % regress trial initiation time vs. previous reward reward
    [betas.latency(rr,:), tau.latency(rr)] =...
        regress_latency_vs_rew(A, nback, false, true);
end

% Compare conditional wait/trial initaition time
pPrevVolWT = signrank(wtByPrevRew(:,1), wtByPrevRew(:,2));
pPrevVolTI = signrank(latByPrevRew(:,1), latByPrevRew(:,2));

% Compare taus fit to previous reward regression coefficients for wait
% time and trial initation time
pTITau = signrank(tau.wt, tau.latency);

% Correlate taus fit to previous reward regression coefficients for wait
% time and trial initation time
[r, pTauCorr] = corr(tau.latency, tau.wt, rows='complete');

% Pull fast and slow taus for latency
first_q = tau.latency < prctile(tau.latency, 100-N);
last_q = tau.latency > prctile(tau.latency, N);

% Perform shuffle test- comparing regression coefficients for both wait
% time and trial initiation time for fast and slow latency rats
[~, ~, pvalLat] =...
    myTauShuffleTest(first_q, last_q, betas.latency, Nshuffle);
[~, ~, pvalWT] =...
    myTauShuffleTest(first_q, last_q, betas.wt, Nshuffle);

%% Plot main figure
% Height (ht) and width (wd) of figure
ht = 3;
wd = 4;

clear l
figure

%--------------------------------------------------------------------------
% 2A - WT dynamics
%--------------------------------------------------------------------------

subplot(ht, wd, 1); hold on
l(7) = shadedErrorBar(-twin:twin, mean(WT.ltom, 'omitnan'),...
    std(WT.ltom, 'omitnan')./sqrt(length(ratList)),...
    LineProps={'color', 'b'});
l(8) = shadedErrorBar(-twin:twin, mean(WT.htom, 'omitnan'),...
    std(WT.htom, 'omitnan')./sqrt(length(ratList)),...
    LineProps={'color', 'r'});

yl = ylim;
fill([0 twin twin 0], [yl(1) yl(1) yl(2) yl(2)],...
    'k', facealpha=0.1, linestyle='none')

ylim(yl)

xlabel('Trials from block switch')
title({'Wait-time', 'Adapt -> Mixed'})

%--------------------------------------------------------------------------
% 2B - TI Time dynamics
%--------------------------------------------------------------------------

subplot(ht, wd, 2); hold on
l(9) = shadedErrorBar(-twin:twin, mean(Lat.ltom, 'omitnan'),...
    std(Lat.ltom, 'omitnan')./sqrt(length(ratList)),...
    LineProps={'color', 'b'});
l(10) = shadedErrorBar(-twin:twin, mean(Lat.htom, 'omitnan'),...
    std(Lat.htom, 'omitnan')./sqrt(length(ratList)),...
    LineProps={'color', 'r'});

yl = ylim;
fill([0 twin twin 0], [yl(1) yl(1) yl(2) yl(2)],...
    'k', facealpha=0.1, linestyle='none')

title({'Trial Initiation Time', 'Adapt -> Mixed'})

%--------------------------------------------------------------------------
% Trial initiation time coeff
%--------------------------------------------------------------------------

subplot(ht, wd, 3); hold on
plot(1:6, betas.latency(:, 2:7),...
    linewidth=0.5, color=[0.8 0.8 0.8])
l(13) = shadedErrorBar(1:6, mean(betas.latency(:, 2:7)),...
    std(betas.latency(:, 2:7))./sqrt(size(betas.latency, 1)),...
    lineprops={'color', lat_color, 'linewidth', 3});

% yline(0, 'k--', linewidth=1, alpha=1)
xticks(0:6)

% ylim([-0.07, 0.015])

xlabel('Trials Back')
ylabel('Regression Coefficient')

%--------------------------------------------------------------------------
% Wait time coeff
%--------------------------------------------------------------------------

subplot(ht, wd, 4); hold on
plot(0:6, betas.wt(:, 1:7),...
    linewidth=0.5, color=[0.8 0.8 0.8])
l(14) = shadedErrorBar(0:6, mean(betas.wt(:, 1:7)),...
    std(betas.wt(:, 1:7))./sqrt(size(betas.wt, 1)),...
    lineprops={'color', wt_color, 'linewidth', 3});

% yline(0, 'k--', linewidth=1, alpha=1)

xticks(0:6)

xlim([-0.5, 6.5])
ylim([-0.13, 0.4])

xlabel('Trials Back')
ylabel('Regression Coefficient')

%--------------------------------------------------------------------------
% Tau histo
%--------------------------------------------------------------------------

subplot(ht, wd, 5); hold on
h = histogram(tau.latency, FaceColor=lat_color,...
    binedges=linspace(0, 5, 25));
histogram(tau.wt, FaceColor=wt_color, binedges=h.BinEdges);

xline(mean(tau.wt, 'omitnan'),...
    '--', color=wt_color, linewidth=2, alpha=1)
xline(mean(tau.latency, 'omitnan'),...
    '--', color=lat_color, linewidth=2, alpha=1)

xlim([-0.15 5.15])
title(pTITau)

%--------------------------------------------------------------------------
% WT vs. TI time taus
%--------------------------------------------------------------------------

subplot(ht, wd, 6); hold on
plot(tau.latency, tau.wt, 'k.', markersize=10)
l1 = refline(1, 0);
set(l1, LineStyle='--', Color=[0 0 0], LineWidth=1)

xlabel('\tau latency')
ylabel('\tau wait time')

xlim([-0.25 5.5])
ylim([-0.25 5.5])

xlabel('Trial Initiation Time \tau')
ylabel('Wait Time \tau')

title(['r = ' num2str(r) ', p = ' num2str(pTauCorr)])

%--------------------------------------------------------------------------
% Latency Tau Percentiles - Latency
%--------------------------------------------------------------------------

subplot(ht, wd, 7); hold on
plot(1:9,  betas.latency(first_q, 2:end-1),...
    linewidth=0.25, color=[firstQColor 0.15])
plot(1:9,  betas.latency(last_q, 2:end-1),...
    linewidth=0.25, color=[lastQColor 0.15])

l(1) = shadedErrorBar(1:9, mean(betas.latency(first_q, 2:end-1)),...
    std(betas.latency(first_q,2:end-1))./sqrt(sum(first_q)),...
    'lineProps', {'-', 'color', firstQColor});
l(2) = shadedErrorBar(1:9, mean(betas.latency(last_q, 2:end-1)),...
    std(betas.latency(last_q,2:end-1))./sqrt(sum(last_q)),...
    'lineProps', {'-', 'color', lastQColor});
% yline(0, 'k--', linewidth=1, alpha=1)

l(1).mainLine.DisplayName = 'Fast Tau';
l(2).mainLine.DisplayName = 'Slow Tau';

xlim([0 10])
% ylim([-0.04 0.005])
xticks(1:2:9)

legend([l(1).mainLine, l(2).mainLine], location='best')

title(['Trial Initiation Time, p = ' num2str(pvalLat)])

%--------------------------------------------------------------------------
% Latency Tau Percentiles - Wait time
%--------------------------------------------------------------------------

subplot(ht, wd, 8); hold on
plot(0:9,  betas.wt(first_q, 1:end-1),...
    linewidth=0.25, color=[firstQColor 0.15])
plot(0:9,  betas.wt(last_q, 1:end-1),...
    linewidth=0.25, color=[lastQColor 0.15])

l(3) = shadedErrorBar(0:9, mean(betas.wt(first_q,1:end-1)),...
    std(betas.wt(first_q,1:end-1))./sqrt(sum(first_q)),...
    'lineProps', {'-', 'color', firstQColor});
l(4) = shadedErrorBar(0:9, mean(betas.wt(last_q,1:end-1)),...
    std(betas.wt(last_q,1:end-1))./sqrt(sum(last_q)),...
    'lineProps', {'-', 'color', lastQColor});
% yline(0, 'k--', linewidth=1, alpha=1)

xlim([-0.5 7])
xticks(1:2:7)

title(['Wait Time, p = ' num2str(pvalWT)])

l(3).mainLine.DisplayName = 'Fast Tau';
l(4).mainLine.DisplayName = 'Slow Tau';

arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), l)
set(gcf, Position=[555 284 1126 663])

%--------------------------------------------------------------------------
% Model agnostic cartoon
%--------------------------------------------------------------------------

subplot(ht, wd, 9); hold on
plot(1:2, [1 1], 'o-',...
    color=[0.3 0.3 0.3], MarkerFaceColor=[0.3 0.3 0.3]) % mb
plot(1:2, [1.5 0.5], 'o-',...
    color=[0.7 0.7 0.7], MarkerFaceColor=[0.7 0.7 0.7]) % mf

ylabel('\Delta Wait Time')

legend('Model-Based', 'Model-Free')
xlim([0.5 2.5])
ylim([0 2])

xticks(1:2)
xticklabels({'< 20', '> 20'})
xlabel('Previous Offer')

%--------------------------------------------------------------------------
% wt by prev rew
%--------------------------------------------------------------------------

subplot(ht, wd, 10); hold on
arrayfun(@(rr)...
    plot([1 2], [wtByPrevRew(rr,1) wtByPrevRew(rr,2)],...
    linewidth=0.5, color=[0.8 0.8 0.8]),...
    1:length(ratList))
plot([1 2], mean(wtByPrevRew),...
    linewidth=2, color=mbColor)

xlim([0.5 2.5])
ylim([-1, 2])

xticks(1:2)
xticklabels({'< 20', '> 20'})
xlabel('Previous Offer')
ylabel('Wait time (z-score)')

title(['p = ' num2str(pPrevVolWT)])

%--------------------------------------------------------------------------
% iti by prev rew
%--------------------------------------------------------------------------

subplot(ht, wd, 11); hold on
arrayfun(@(rr)...
    plot([1 2], [latByPrevRew(rr,1) latByPrevRew(rr,2)],...
    linewidth=0.5, color=[0.8 0.8 0.8]),...
    1:length(ratList))
plot([1 2], mean(latByPrevRew),...
    linewidth=2, color=mfColor)

title(['p = ' num2str(pPrevVolTI)])
xlim([0.75 2.25])
ylim([-0.3, 0.2])

xticks(1:2)
xticklabels({'< 20', '> 20'})
xlabel('Previous Offer')
ylabel('Trial initiation time (z-score)')

%--------------------------------------------------------------------------
% hist - rats
%--------------------------------------------------------------------------

wt_diff = wtByPrevRew(:,1)-wtByPrevRew(:,2);
lat_diff = latByPrevRew(:,1)-latByPrevRew(:,2);

binwidth = 0.1;
subplot(ht, wd, 12); hold on
h1 = histogram(lat_diff,...
    EdgeColor=mfColor, displaystyle='stair',...
    LineWidth=2, binWidth=binwidth, BinLimits=[-1 1]);
histogram(lat_diff(pLat<0.05),...
    FaceColor=mfColor, LineStyle='none', facealpha=1,...
    BinEdges=h1.BinEdges);

h2 = histogram(wt_diff,...
    EdgeColor=mbColor, displaystyle='stair',...
    LineWidth=2, binWidth=binwidth, BinLimits=[-1 1]);
histogram(wt_diff(pWt<0.05),...
    FaceColor=mbColor, LineStyle='none', facealpha=1,...
    BinEdges=h2.BinEdges);

xline(0, 'k-', LineWidth=2, Alpha=1)
xline(mean(wt_diff), '--', color=mbColor, LineWidth=2, Alpha=1)
xline(mean(lat_diff), '--', color=mfColor, LineWidth=2, Alpha=1)

title(signrank(wt_diff, lat_diff))

set(gcf, Position=[596 81 1085 866], renderer='painters');
arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), l)


%--------------------------------------------------------------------------
% Insets - cartoon
%--------------------------------------------------------------------------
figure
xs = -2:0.01:2;
subplot(ht, wd, 1); hold on
plot(xs, normpdf(xs, 0, 0.1), color=[0.3 0.3 0.3])
plot(xs, normpdf(xs, 1, 0.1), color=[0.7 0.7 0.7])
xline(0, 'r--', linewidth=2, alpha=1)

xlim([-2 2])
ylim([0.05, 4])

ylabel('rats')
yticks([])
xticks(0)
set(gcf, Position=[596 81 1085 866], renderer='painters');

subplot(ht, wd, 2)
shadedErrorBar(1:7, mean(betas.wt(first_q,2:8)),...
    std(betas.wt(first_q,2:8))./sqrt(sum(first_q)),...
    'lineProps', {'-', 'color', firstQColor});
shadedErrorBar(1:7, mean(betas.wt(last_q,2:8)),...
    std(betas.wt(last_q,2:8))./sqrt(sum(last_q)),...
    'lineProps', {'-', 'color', lastQColor});
xticks(1:2:7)


end