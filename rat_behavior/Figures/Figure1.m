function [WT, Lat] = Figure1(datadir, codedir)
%Figure1 - Wait time and trial initiation time were modulated by the value 
% of the environment
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

% example rat names for data
examprat_wt_curve = 'S038';
examprat_wt_dist = 'J040';
examprat_lat = 'J031';

% color scheme for figure
wt_color = '#a16ae8';
iti_color = '#50C878';

% Preallocate matricies
WT = struct;
[WT.high, WT.low, WT.mixed]  = deal(nan(length(ratList), 5));

Lat = struct;
[Lat.high, Lat.low,...
    ps, ntrials]=...
    deal(nan(length(ratList), 2));

% Process each rat's data
for rr = 1:length(ratList)
    fprintf('%d out of %d\n', rr, length(ratList));

    % Load Data
    fname = strcat(['ratTrial_', ratList{rr}, '.mat']);
    a = load([datadir 'A_Structs_Final' filesep fname]);
    A = a.A;

    % Average wait time and SEM as a function of reward in each block
    [high, low, mix, ps(rr,1)] = wtcurves(A);

    WT.high(rr,:) = high.wt;
    WT.low(rr,:) = low.wt;
    WT.mixed(rr,:) = mix.wt;

    if strcmp(ratList{rr}, examprat_wt_curve) % for example rat
        WT.example_rat.mean  = [mix.wt; high.wt; low.wt];
        WT.example_rat.error = [mix.er; high.er; low.er];
    elseif strcmp(ratList{rr}, examprat_wt_dist)
        WT.example_rat.raw = A.wait_time(A.catch & A.optout);
    end

    % Average trial initation time and SEM in each block
    [hi, lo, ~, ps(rr,2)] = iticurves(A);
    
    Lat.raw.high(rr,:) = hi.raw(1);
    Lat.raw.low(rr,:) = lo.raw(1);

    Lat.high(rr,:) = hi.z(1);
    Lat.low(rr,:) = lo.z(1);

    if strcmp(ratList{rr}, examprat_wt_curve) % for example rat
        Lat.example_rat  = [lo.raw; hi.raw];
    end

    % Number of sessions and trials
    ntrials(rr,:) = [length(A.date) length(A.reward)];
end

% Compare wait time for 20 uL in high vs. low blocks
pwt = signrank(WT.high(:,3), WT.low(:,3));
wt_ratio = WT.high(:,3)./WT.low(:,3);

% Compare trial initiation time in high vs. low blocks
piti = signrank(Lat.raw.high(:,1), Lat.raw.low(:,1));
lat_ratio = Lat.raw.high(:,1)./Lat.raw.low(:,1);

% calculate binwidths for figure
[~, edges_iti] = histcounts(lat_ratio, NumBins=15);
iti_binwidth = mean(diff(edges_iti));
wt_binwidth = (0.27/0.8)*iti_binwidth;

fprintf(['Median sessions: %d\n'...
    'Median trials: %d\n'],...
    median(ntrials(:,1)), median(ntrials(:,2)))

%% Plot main figure
clear l

wt = 4;
ht = 2;

alpha = 0.05;

figure
%--------------------------------------------------------------------------
% Example rat dist
%--------------------------------------------------------------------------
subplot(ht, wt, 1)
histogram(WT.example_rat.raw, DisplayStyle='stair',...
    edgecolor='k', linewidth=1)
xline(mean(WT.example_rat.raw, 'omitnan'), 'k', alpha=1, linewidth=1)

xlim([-2 31])
ylim([0 630])

xlabel('Wait time (s)')
ylabel('N (trials)')

title(['Rat ' examprat_wt_dist])
set(gca, Box='off')

%--------------------------------------------------------------------------
% Block fig
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Example rat curves
%--------------------------------------------------------------------------

xs = 1:5;
subplot(ht, wt, 3)
l(1) = shadedErrorBar(1:5, WT.example_rat.mean(1,:),...
    WT.example_rat.error(1,:), LineProps={'k'});
l(2) = shadedErrorBar(1:5, WT.example_rat.mean(2,:),...
    WT.example_rat.error(2,:), LineProps={'r'});
l(3) = shadedErrorBar(1:5, WT.example_rat.mean(3,:),...
    WT.example_rat.error(3,:), LineProps={'b'});

xlim([0.5 5.5])
ylim([8.15 13.5])

xticks(1:5)
xticklabels({'5', '10', '20', '40', '80'})

xlabel('Offered reward')
ylabel('Mean wait time (s)')
title(['Rat ' examprat_wt_curve])

%--------------------------------------------------------------------------
% Population wt curves
%--------------------------------------------------------------------------

subplot(ht, wt, 4); hold on
l(4) = shadedErrorBar(xs, mean(WT.mixed),...
    std(WT.mixed)./sqrt(length(ratList)));

l(5) = shadedErrorBar(xs, mean(WT.high),...
    std(WT.high)./sqrt(length(ratList)),...
    LineProps={'r'});

l(6) = shadedErrorBar(xs, mean(WT.low),...
    std(WT.low)./sqrt(length(ratList)),...
    LineProps={'b'});

xlim([0.5 5.5])
ylim([10.25 14.3])

xticks(1:5)
xticklabels({'5', '10', '20', '40', '80'})

xlabel('Offered reward')
ylabel('Mean wait time (s)')
title(['N = ' num2str(length(ratList))])

%--------------------------------------------------------------------------
% WT ratio histogram
%--------------------------------------------------------------------------

subplot(ht, wt, 5); hold on
h1 = histogram(wt_ratio, facecolor='none', BinWidth=wt_binwidth);
histogram(wt_ratio(ps(:,1)<alpha),...
    facecolor=wt_color, BinEdges = h1.BinEdges)

xline(1, 'k--', linewidth=1, alpha=1)
xline(mean(wt_ratio), 'k', linewidth=1, alpha=1)

xlabel('Wait time ratio (20 High/Low)')
ylabel('N (rats)')

xlim(1+[-0.27, 0.27])
ylim([0 42])
title(pwt)

%--------------------------------------------------------------------------
% Example rat block lat
%--------------------------------------------------------------------------

subplot(ht, wt, 6); hold on
errorbar(1, Lat.example_rat(1,1), Lat.example_rat(1,2),...
    'bo', markerfacecolor='b', linewidth=1, capsize=0)

errorbar(2, Lat.example_rat(2, 1), Lat.example_rat(2,2),...
    'ro', markerfacecolor='r', linewidth=1, capsize=0)

xlim([0.5 2.5])
ylim([2, 4.5])

xticks(1:2)
xticklabels({'Low', 'High'})

xlabel('Reward block')
ylabel('Mean trial initiation time (s)')

title(['Rat ' examprat_lat])

%--------------------------------------------------------------------------
% Population block lat
%--------------------------------------------------------------------------

subplot(ht, wt, 7); hold on
arrayfun(@(r) plot([1 2], [Lat.low(r,1) Lat.high(r,1)],...
    color=[0.7 0.7 0.7], LineWidth=0.25), 1:length(ratList))
plot([1 2], [mean(Lat.low(:,1)) mean(Lat.high(:,1))], 'k')

xlim([0.5 2.5])
ylim([-0.23 0.23])

xticks(1:2)
xticklabels({'Low', 'High'})

xlabel('Reward block')
ylabel('Mean trial initiation time (z-score)')

title(['N = ' num2str(length(ratList))])

%--------------------------------------------------------------------------
% Lat ratio histogram
%--------------------------------------------------------------------------

subplot(ht, wt, 8); hold on
h2 = histogram(lat_ratio, facecolor='none', Binwidth=iti_binwidth);
histogram(lat_ratio(ps(:,2)<alpha),...
    facecolor=iti_color, BinEdges = h2.BinEdges)
xline(1, 'k--', linewidth=1, alpha=1)
xline(mean(lat_ratio), 'k', linewidth=1, alpha=1)

xlim(1+[-0.8, 0.8])

xlabel('Trial initiation time ratio')
ylabel('N (rats)')

title(piti)

set(gcf, Position=[162 307 1519 638])
arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), l)
end