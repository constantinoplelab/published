function FigureS12(datadir, codedir)
%FigureS12 Males and females have comparable wait time ratios (top) and 
% trial initiation time ratios (bottom)
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
sexList = a.sexList; % corresponding sex

% block transition dynamic parameters
twin = 20; % number of trials around each block to pull
smoothfactor = 10; % size of smoothing window

% Preallocate matricies
WT = struct;
[WT.high, WT.low, WT.mixed]  = deal(nan(length(ratList), 5));

Lat = struct;
[Lat.high, Lat.low]=...
    deal(nan(length(ratList), 2));

[WT.ltom, WT.htom, Lat.ltom, Lat.htom] =...
    deal(nan(length(ratList), 2*twin+1));

% Process each rat's data
for rr = 1:length(ratList)
    fprintf('%d out of %d\n', rr, length(ratList));

    % Load Data
    fname = strcat(['ratTrial_', ratList{rr}, '.mat']);
    a = load([datadir 'A_Structs_Final' filesep fname]);
    A = a.A;

    % Average wait time and SEM as a function of reward in each block
    [high, low, mix] = wtcurves(A);

    WT.high(rr,:) = high.wt;
    WT.low(rr,:) = low.wt;
    WT.mixed(rr,:) = mix.wt;

    % Average trial initation time and SEM in each block
    [hi, lo] = iticurves(A);

    Lat.high(rr,:) = hi.raw(1);
    Lat.low(rr,:) = lo.raw(1);

    % Wait time dynamics around block transitions
    [WT.ltom(rr,:), WT.htom(rr,:)] =...
        block_dynamics_wt(A, twin, smoothfactor);

    % trial intiation time dynamics around block transitions
    [Lat.ltom(rr,:), Lat.htom(rr,:)] =...
        block_dynamics_latency(A, twin, smoothfactor);
end

wt_ratio = WT.high(:,3)./WT.low(:,3);
lat_ratio = Lat.high(:,1)./Lat.low(:,1);

maleRats = strcmpi(sexList, 'M');
femaleRats = strcmpi(sexList, 'F');

pWT = ranksum(wt_ratio(maleRats), wt_ratio(femaleRats));
pTI = ranksum(lat_ratio(maleRats), lat_ratio(femaleRats));

[nWTMales, edgesWT] = histcounts(wt_ratio(maleRats), 17);
nWTFemales = histcounts(wt_ratio(femaleRats), edgesWT);
xsWT = (edgesWT(1:end-1) + edgesWT(2:end))/2;

[nTIMales, edgesTI] = histcounts(lat_ratio(maleRats), 17);
nTIFemales = histcounts(lat_ratio(femaleRats), edgesTI);
xsTI = (edgesTI(1:end-1) + edgesTI(2:end))/2;

%% Plot figure

figure
subplot(2, 2, 1); hold on
bar(xsWT, nWTMales);
bar(xsWT, -nWTFemales);

xline(mean(wt_ratio(maleRats)), 'b--', linewidth=1, alpha=1)
xline(mean(wt_ratio(femaleRats)), 'r--', linewidth=1, alpha=1)
xline(1, 'k--', linewidth=1, alpha=1)

ylabel('n (rats)')
xlabel('Wait time ratio')

title({['Wait time ' num2str(pWT)], [sum(maleRats) sum(femaleRats)]})

subplot(2, 2, 3); hold on
bar(xsTI, nTIMales);
bar(xsTI, -nTIFemales);
title(pTI)

xline(mean(lat_ratio(maleRats)), 'b--', linewidth=1, alpha=1)
xline(mean(lat_ratio(femaleRats)), 'r--', linewidth=1, alpha=1)
xline(1, 'k--', linewidth=1, alpha=1)

ylabel('n (rats)')
xlabel('Trial initiation time ratio')

title(['Trial initiation time ' num2str(pTI)])


subplot(2, 2, 2); hold on

a = fill([-twin 0 0 -twin], [-0.2 -0.2 0.3 0.3],...
    'k', facealpha=0.1, linestyle='none');

l(1) = shadedErrorBar(-twin:twin,...
    mean(WT.ltom(maleRats,:)),...
    std(WT.ltom(maleRats,:))./sqrt(sum(maleRats)),...
    lineprops={'color', 'b', });
l(2) = shadedErrorBar(-twin:twin,...
    mean(WT.htom(maleRats,:)),...
    std(WT.htom(maleRats,:))./sqrt(sum(maleRats)),...
    lineprops={'color', 'r', });

[l(1).mainLine.DisplayName, l(2).mainLine.DisplayName] = deal('Male');

l(3) = shadedErrorBar(-twin:twin,...
    mean(WT.ltom(femaleRats,:)),...
    std(WT.ltom(femaleRats,:))./sqrt(sum(femaleRats)),...
    lineprops={'--b', });
l(4) = shadedErrorBar(-twin:twin,...
    mean(WT.htom(femaleRats,:)),...
    std(WT.htom(femaleRats,:))./sqrt(sum(femaleRats)),...
    lineprops={'--r', });

[l(3).mainLine.DisplayName, l(4).mainLine.DisplayName] = deal('Female');
legend()

subplot(2, 2, 4); hold on
fill([-twin 0 0 -twin], [-0.15 -0.15 0.2 0.2],...
    'k', facealpha=0.1, linestyle='none')

shadedErrorBar(-twin:twin,...
    mean(Lat.ltom(maleRats,:)),...
    std(Lat.ltom(maleRats,:))./sqrt(sum(maleRats)),...
    lineprops={'color', 'b', })
shadedErrorBar(-twin:twin,...
    mean(Lat.htom(maleRats,:)),...
    std(Lat.htom(maleRats,:))./sqrt(sum(maleRats)),...
    lineprops={'color', 'r', })

shadedErrorBar(-twin:twin,...
    mean(Lat.ltom(femaleRats,:)),...
    std(Lat.ltom(femaleRats,:))./sqrt(sum(femaleRats)),...
    lineprops={'--b', })
shadedErrorBar(-twin:twin,...
    mean(Lat.htom(femaleRats,:)),...
    std(Lat.htom(femaleRats,:))./sqrt(sum(femaleRats)),...
    lineprops={'--r', })


set(gcf, position=[852 80 829 865])
shg

end