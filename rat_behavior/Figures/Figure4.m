function [] = Figure4(datadir, codedir)
%Figure4 - Block sensitivity for wait times requires structure learning
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
ratList = a.ratList; % List of rats

% block transition dynamic parameters
twin = 20; % number of trials around each block to pull
smoothfactor = 10; % size of smoothing window

% Parameter for example sessions for high to mixed transition
nsess = 10; % Pull first nsess and last nsess to average across rats
chunkSize = 5; % Number of sessions to combine

% Preallcate matricies
[htom_iti_early, ltom_iti_early, mtol_iti_early, mtoh_iti_early,...
    htom_iti_late, ltom_iti_late, mtol_iti_late, mtoh_iti_late,...
    htom_wt_early, ltom_wt_early, mtol_wt_early, mtoh_wt_early,...
    htom_wt_late, ltom_wt_late, mtol_wt_late, mtoh_wt_late] =...
    deal(nan(length(ratList), 2*twin+1));

[data_early, data_late] = deal(cell(1, length(ratList)));

[ltomChunk, htomChunk,...
    mtolChunk, mtohChunk] = deal(cell(1, length(ratList)));

% Loop over rats
for rr = 1:length(ratList)    
    % Load S Struct- has data over training that gets removed in A struct
    s = load([datadir filesep 'S_Structs' filesep...
        'SStruct_' ratList{rr} '.mat']);
    S_all = s.S;

    % Pull the training stage
    trainingStage = cellfun(@(sess) sess.TrainingStage(1), S_all.pd);

    % Pull the protocol - very early rats ran on a slightly different
    % protocol so we don't want to include them (normally gets filtered but
    % this is raw data)
    protocol = cellfun(@(sess) sess.Protocol, S_all.pd,...
        UniformOutput=false);

    % Pull final stage (paper says there are 8 stages, but final stage in
    % data is coded as 9 because of a holdover from a previous iteration of
    % the task that's not used). Also remove experimental rats
    usethese = trainingStage==9 &...
        contains(protocol, 'RWTautowait') &...
        ~contains(protocol, 'Opto');
    S = struct('pd', {S_all.pd(usethese)},...
        'peh', {S_all.peh(usethese)});

    % If the rat has no data, skip it
    if isempty(S.pd)
        continue
    end

    % Pull early vs. late sessions
    S_early =...
        struct('pd', {S.pd(1:nsess)},...
        'peh', {S.peh(1:nsess)});
    data_early{rr} = parse_data_from_mysql(S_early, [], 'all');

    S_late =...
        struct('pd', {S.pd(end-nsess+1:end)},...
        'peh', {S.peh(end-nsess+1:end)});
    data_late{rr} = parse_data_from_mysql(S_late, [], 'all');

    % Pull trial initiation time dynamics for early vs. late
    [ltom_iti_early(rr,:), htom_iti_early(rr,:),...
        mtol_iti_early(rr,:), mtoh_iti_early(rr,:)] =...
        block_dynamics_latency(data_early{rr}, twin, smoothfactor);
    [ltom_iti_late(rr,:), htom_iti_late(rr,:),...
        mtol_iti_late(rr,:), mtoh_iti_late(rr,:)] =...
        block_dynamics_latency(data_late{rr}, twin, smoothfactor);

    [ltom_wt_early(rr,:), htom_wt_early(rr,:),...
        mtol_wt_early(rr,:), mtoh_wt_early(rr,:)] =...
        block_dynamics_wt(data_early{rr}, twin, smoothfactor);
    [ltom_wt_late(rr,:), htom_wt_late(rr,:),...
        mtol_wt_late(rr,:), mtoh_wt_late(rr,:)] =...
        block_dynamics_wt(data_late{rr}, twin, smoothfactor);

    % Look at dynamics (specificaly overshoot) across sessions. Chunk
    % sessions to have enough data
    numSessions = length(S.pd);
    numChunk = ceil(numSessions/5);
    Is = reshape([1:numSessions,...
        nan(1, mod(chunkSize-mod(numSessions, chunkSize), chunkSize))],...
        [chunkSize, numChunk]);

    [ltomChunk{rr}, htomChunk{rr},...
        mtolChunk{rr}, mtohChunk{rr}] = deal(nan(numChunk, 2*twin+1));

    for cc = 1:numChunk
        fprintf('%d out of %d: %.2f\n',...
            rr, length(ratList), cc/numChunk)

        usethese = Is(~isnan(Is(:,cc)), cc);

        SChunk = struct('pd', {S.pd(usethese)},...
            'peh', {S.peh(usethese)});
        dataChunk = parse_data_from_mysql(SChunk, [], 'all');

        [ltomChunk{rr}(cc,:), htomChunk{rr}(cc,:),...
            mtolChunk{rr}(cc,:), mtohChunk{rr}(cc,:)] =...
            block_dynamics_latency(dataChunk, twin, smoothfactor);
    end
end
%%
max_sess = max(cellfun(@length, ltomChunk));
pool_htom = cell(1, max_sess);
htomAvg = nan(max_sess, 2*twin+1);

for ss = 1:max_sess
    pool_htom{ss} = nan(length(ratList), 2*twin+1);

    for rr = 1:length(ratList)
        disp([ss rr])
        if ss > size(htomChunk{rr}, 1)
            pool_htom{ss}(rr,:) = nan(1, 2*twin+1);
        else
            pool_htom{ss}(rr,:) = htomChunk{rr}(ss,:);
        end

    end

    htomAvg(ss,:) = mean(pool_htom{ss}, 'omitnan');
end

overshoot = max(htomAvg(:, twin+1:end), [], 2) - htomAvg(:, end);

%%
nsess = 50;

[B1, B2,...
    B1iti, B2iti,...
    dwt1, dwt2,...
    diti1, diti2] =...
    macro_regress_overtraining_V2(nsess, ratList,...
    [datadir 'A_structs_all' filesep]);

%%
clear l

grp = 2;
nsess = 50;

xvals = grp:grp:nsess;
lprops = {'k'; 'b'; 'r'};

ht = 4;
wd = 2;

figure
%--------------------------------------------------------------------------
% Wait Time
%--------------------------------------------------------------------------

subplot(ht, wd, 1); hold on
shadedErrorBar(xvals,...
    mean(dwt1, 'omitnan'),...
    std(dwt1, 'omitnan')./sqrt(sum(~isnan(dwt1))),...
    'lineprops', {'color', '#8B6DB0', 'linestyle', 'none'})

shadedErrorBar(xvals,...
    mean(diti1, 'omitnan'),...
    std(diti1, 'omitnan')./sqrt(sum(~isnan(diti1))),...
    'lineprops', {'color', '#5DBD77', 'linestyle', 'none'})

yline(1, 'k', linewidth=2, alpha=1)

xlim([0 nsess]);
ylim([0.55 1.1]);

xlabel('Training session');
ylabel('Wait time ratio');
title(strcat(['First ', num2str(nsess), ' sessions']));

subplot(ht, wd, 2); hold on
shadedErrorBar(xvals,...
    mean(dwt2, 'omitnan'),...
    std(dwt2, 'omitnan')./sqrt(sum(~isnan(dwt2))),...
    'lineprops', {'color', '#8B6DB0', 'linestyle', 'none'})
shadedErrorBar(xvals,...
    mean(diti2, 'omitnan'),...
    std(diti2, 'omitnan')./sqrt(sum(~isnan(diti2))),...
    'lineprops', {'color', '#5DBD77', 'linestyle', 'none'})

yline(1, 'k', linewidth=2, alpha=1)

xlim([0 nsess]);
ylim([0.55 1.1]);

xlabel('Training session');
ylabel('Wait time ratio');
title(strcat(['Last ', num2str(nsess), ' sessions']));

myc = [0 0 0; 0 0 1; 1 0 0];
for k = 1:3
    if k == 3
        i=-1;
    else
        i=1;
    end
    
    subplot(ht, wd, 3); hold on
    shadedErrorBar(xvals,...
        i*mean(B1{k}, 'omitnan'),...
        std(B1{k}, 'omitnan')./sqrt(sum(~isnan(B1{k}))),...
        'lineprops',...
        {'color', lprops{k}, 'linestyle', 'none'}); %#ok<*SAGROW>
    
    [~, yfit] = get_tau(xvals, i*mean(B1{k}, 'omitnan'));
    plot(xvals, yfit, 'Color', myc(k,:), 'LineStyle', '--');
    line([0 nsess], [0 0], 'Color', [0 0 0]);
    
    xlim([0 nsess]);
    ylim([-0.05 .3]);
    xlabel('Training session');
    ylabel('Regression coefficient');
    
    subplot(ht, wd, 4); hold on
    shadedErrorBar(xvals,...
        i*mean(B2{k}, 'omitnan'),...
        std(B2{k}, 'omitnan')./sqrt(sum(~isnan(B2{k}))),...
        'lineprops', {'color', lprops{k}, 'linestyle', 'none'});
    
    [~, yfit] = get_tau(xvals, i*mean(B2{k}, 'omitnan'));
    plot(xvals, yfit, 'Color', myc(k,:), 'LineStyle', '--');
    
    line([0 nsess], [0 0], 'Color', [0 0 0]);
    
    xlim([0 nsess]);
    ylim([-0.05 .3]);
    xlabel('Training session');
    ylabel('Regression coefficient');
end

for k = 1:3
    subplot(ht, wd, 5); hold on
    shadedErrorBar(xvals,...
        mean(B1iti{k}, 'omitnan'),...
        std(B1iti{k}, 'omitnan')./sqrt(sum(~isnan(B1iti{k}))),...
        'lineprops', {'color', lprops{k}, 'linestyle', 'none'})
    
    [~, yfit] = get_tau(xvals, mean(B1iti{k}, 'omitnan'));
    plot(xvals, yfit, 'Color', myc(k,:), 'LineStyle', '--');
    line([0 nsess], [0 0], 'Color', [0 0 0]);
    
    xlim([0 nsess]);
    ylim([-8e-3 1e-4]);
    xlabel('Training session');
    ylabel('Regression coefficient');
    
    subplot(ht, wd, 6); hold on
    shadedErrorBar(xvals,...
        mean(B2iti{k}, 'omitnan'),...
        std(B2iti{k}, 'omitnan')./sqrt(sum(~isnan(B2iti{k}))),...
        'lineprops', {'color', lprops{k}, 'linestyle', 'none'})
    
    [~, yfit] = get_tau(xvals, mean(B2iti{k}, 'omitnan'));
    plot(xvals, yfit, 'Color', myc(k,:), 'LineStyle', '--');
    
    line([0 nsess], [0 0], 'Color', [0 0 0]);
    
    xlim([0 nsess]);
    ylim([-8e-3 1e-4]);
    xlabel('Training session');
    ylabel('Regression coefficient');
end

%--------------------------------------------------------------------------
% Overshoot
%--------------------------------------------------------------------------
subplot(ht, wd, 7)
shadedErrorBar(-twin:twin,...
    mean(htom_iti_early),...
    std(htom_iti_early)./sqrt(size(htom_iti_early, 1)),...
    lineprops={'color', 'r'});
shadedErrorBar(-twin:twin,...
    mean(htom_iti_late, 'omitnan'),...
    std(htom_iti_late, 'omitnan')./sqrt(size(htom_iti_late, 1)),...
    lineprops={'color', 'r', 'linestyle', '--'});

subplot(ht, wd, 8)
plot(5:5:max_sess*5, smooth(overshoot, 10))

xlim([0 60])
set(gca, Box='off')

xlabel('session')
ylabel('overshoot')

set(gcf, Position=[1194 80 487 865])
end