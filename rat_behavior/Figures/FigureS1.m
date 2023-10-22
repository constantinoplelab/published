function FigureS1(datadir, codedir)
%FigureS1 -  Trial initiation time in units of seconds 
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

%% Analyze data

a = load([datadir filesep 'ratList.mat']);
ratList = a.ratList; % list of rats

% Preallocate matricies
RawLat = nan(length(ratList), 2);

% Process each rat's data
for rr = 1:length(ratList)
    fprintf('%d out of %d\n', rr, length(ratList));

    % Load Data
    fname = strcat(['ratTrial_', ratList{rr}, '.mat']);
    a = load([datadir 'A_Structs_Final' filesep fname]);
    A = a.A;

    % Average trial initation time and SEM in seconds for each block
    [hi, lo] = iticurves(A);
    RawLat(rr,:) = [lo.raw(1), hi.raw(1)];
end

%% Plot
pval = signrank(RawLat(:,1), RawLat(:,2));

figure
subplot(1, 2, 1); hold on
plot(1:2, [RawLat(:,1), RawLat(:,2)],...
color=[0.7 0.7 0.7], linewidth=0.25)
plot(1:2, mean([RawLat(:,1), RawLat(:,2)]),...
'k', linewidth=1)

xticks(1:2)
xticklabels({'Low', 'High'})
xlim([0.5 2.5])

ylabel('Trial initiation time (s)')
title(pval)

subplot(1, 2, 2)
histogram(RawLat(:,2) - RawLat(:,1))
xline(mean(RawLat(:,2) - RawLat(:,1)), 'r--', linewidth=1, alpha=1)
xline(0, 'k--', linewidth=1, alpha=1)

ylabel('N (rats)')
xlabel('High - Low')
set(gcf, Position=[980 540 700 420])

shg
end