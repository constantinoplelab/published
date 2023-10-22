function FigureS14(datadir, codedir)
%FigureS14 Wait time curves without threshold have qualitatively similar
% context effects, but longer average wait times.
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

% Preallocate matricies
[mi_wt1, hi_wt1, lo_wt1,...
    mi_wt2, hi_wt2, lo_wt2] = deal(nan(length(ratList), 5));

for rr = 1:length(ratList)
    fprintf('%d out of %d\n', rr, length(ratList));

    % Load Data
    fname = strcat(['ratTrial_', ratList{rr}, '.mat']);
    a = load([datadir 'A_Structs_Final' filesep, fname]);
    A = a.A;

    if isfield(A, 'wait_time_unthresholded')
        usethese = A.optout & A.catch; % Find catch trials

        % Pull Behavioral Data
        wt1 = A.wait_time_unthresholded(usethese); % unthresholded wt
        wt2 = A.wait_time(usethese); % thresholded wt

        rew = convertreward(A.reward(usethese)); % rewards (as 1:5)
        blk = A.block(usethese); % block

        mi_wt1(rr,:) =...
            arrayfun(@(r) mean(wt1(rew==r & blk==1), 'omitnan'), 1:5);
        mi_wt2(rr,:) =...
            arrayfun(@(r) mean(wt2(rew==r & blk==1), 'omitnan'), 1:5);

        hi_wt1(rr,:) =...
            arrayfun(@(r) mean(wt1(rew==r & blk==2), 'omitnan'), 1:5);
        hi_wt2(rr,:) =...
            arrayfun(@(r) mean(wt2(rew==r & blk==2), 'omitnan'), 1:5);

        lo_wt1(rr,:) =...
            arrayfun(@(r) mean(wt1(rew==r & blk==3), 'omitnan'), 1:5);
        lo_wt2(rr,:) =...
            arrayfun(@(r) mean(wt2(rew==r & blk==3), 'omitnan'), 1:5);
    else
        continue
    end
end

%% Plot

figure
subplot(1, 2, 1); hold on
l(1) = shadedErrorBar(1:5, mean(mi_wt2, 'omitnan'),...
    std(mi_wt2, 'omitnan')./sqrt(length(ratList)),...
    lineprops={'color', 'k', 'linewidth', 2});
l(2) = shadedErrorBar(1:5, mean(hi_wt2, 'omitnan'),...
    std(hi_wt2, 'omitnan')./sqrt(length(ratList)),...
    lineprops={'color', 'r', 'linewidth', 2});
l(3) = shadedErrorBar(1:5, mean(lo_wt2, 'omitnan'),...
    std(lo_wt2, 'omitnan')./sqrt(length(ratList)),...
    lineprops={'color', 'b', 'linewidth', 2});

ylabel('Wait time (s)')
xlabel('Reward offer (L)')
title('Treshold')

xlim([0.5 5.5])
ylim([10.25 15.15])

subplot(1, 2, 2); hold on
l(4) = shadedErrorBar(1:5, mean(mi_wt1, 'omitnan'),...
    std(mi_wt1, 'omitnan')./sqrt(length(ratList)),...
    lineprops={'color', 'k', 'linewidth', 2});
l(5) = shadedErrorBar(1:5, mean(hi_wt1, 'omitnan'),...
    std(hi_wt1, 'omitnan')./sqrt(length(ratList)),...
    lineprops={'color', 'r', 'linewidth', 2});
l(6) = shadedErrorBar(1:5, mean(lo_wt1, 'omitnan'),...
    std(lo_wt1, 'omitnan')./sqrt(length(ratList)),...
    lineprops={'color', 'b', 'linewidth', 2});
title('Unthresholded')

xlabel('Reward offer (L)')

xlim([0.5 5.5])
ylim([10.25 15.15])

arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), l)
end