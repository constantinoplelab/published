function FigureS5(datadir, codedir)
%FigureS5 -  Dynamics of wait times and trial initiation times at 
% transitions from mixed to high or low blocks
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
twinLong = 50; % number of trials around each block to pull (longer)
smoothfactor = 10; % size of smoothing window

% Pre-allocate matricies
[WTmtol, WTmtoh, Latmtol, Latmtoh] =...
    deal(nan(length(ratList), 2*twin+1));

[ltomLonger, htomLonger] = deal(nan(length(ratList), 2*twinLong+1));

for rr = 1:length(ratList)
    fprintf('%d out of %d\n', rr, length(ratList));

    % Load Data
    fname = strcat(['ratTrial_', ratList{rr}, '.mat']);
    A = load([datadir 'A_Structs_Final' filesep fname]);
    A = A.A;

% wait time
    % mixed -> high or low dynamics
    [~, ~, WTmtol(rr,:), WTmtoh(rr,:)]=...
        block_dynamics_wt(A, twin, smoothfactor);

% trial initiation time
    % mixed -> high or low dynamics
    [~, ~, Latmtol(rr,:), Latmtoh(rr,:)]=...
        block_dynamics_latency(A, twin, smoothfactor);    

    % larger window
    [ltomLonger(rr,:), htomLonger(rr,:)] =...
        block_dynamics_latency(A, twinLong, smoothfactor);
end

%% Supp. fig - mixed to adapt
 
clear l

figure
subplot(3, 1, 1); hold on

l(1) = shadedErrorBar(-twinLong:twinLong,...
    mean(ltomLonger),...
    std(ltomLonger)./sqrt(length(ratList)), lineprops={'b'});
l(2) = shadedErrorBar(-twinLong:twinLong,...
    mean(htomLonger),...
    std(htomLonger)./sqrt(length(ratList)), lineprops={'r'});

yl = ylim;
fill([0 50 50 0], [yl(1) yl(1) yl(2) yl(2)],...
    'k', facecolor='k', facealpha=0.1, linestyle='none')

ylim(yl)
xlim([-40 40])

ylabel('Delta-Trial initiation time (z-score)')

subplot(3, 1, 2); hold on
l(3) = shadedErrorBar(-twin:twin, mean(WTmtol, 'omitnan'),...
    std(WTmtol, 'omitnan')./sqrt(length(ratList)),...
    LineProps={'color', 'b'});
l(4) = shadedErrorBar(-twin:twin, mean(WTmtoh, 'omitnan'),...
    std(WTmtoh, 'omitnan')./sqrt(length(ratList)),...
    LineProps={'color', 'r'});

yl = ylim;
fill([-twin 0 0 -twin], [yl(1) yl(1) yl(2) yl(2)],...
    'k', facealpha=0.1, linestyle='none')
ylim(yl)

ylabel('Wait time (z-score)')
title({'Wait time', 'Mixed -> Adapt'})

subplot(3, 1, 3); hold on
l(5) = shadedErrorBar(-twin:twin, mean(Latmtol, 'omitnan'),...
    std(Latmtol, 'omitnan')./sqrt(length(ratList)),...
    LineProps={'color', 'b'});
l(6) = shadedErrorBar(-twin:twin, mean(Latmtoh, 'omitnan'),...
    std(Latmtoh, 'omitnan')./sqrt(length(ratList)),...
    LineProps={'color', 'r'});

yl = ylim;
fill([-twin 0 0 -twin], [yl(1) yl(1) yl(2) yl(2)],...
    'k', facealpha=0.1, linestyle='none')
ylim(yl)

xlabel('Trials from block switch')
ylabel('Delta-Trial initiation time (z-score)')

title('Trial initiation time')

set(gcf, Position=[1249 109 432 836])
arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), l)

end