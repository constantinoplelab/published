function [] = Figure5(datadir, codedir)
%Figure5 - Value computations for motivation do not depend on state 
% uncertainty
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

%%

nback = 10;
twin = 20;
smoothfactor = 10;

ratListAll = {'M041', 'M042', 'M043', 'M044', 'M045', 'M046', 'M047',...
    'M048', 'G070', 'E005', 'E008', 'G078', 'G076','G069', 'S042', 'S050'};
startDates = datetime({'2022-10-04', '2022-10-04', '2022-10-04',...
    '2022-10-04', '2022-12-13', '2022-12-13', '2022-12-13',...
    '2022-12-13', '2022-10-24', '2022-10-24', '2022-10-24',...
    '2022-11-07', '2022-11-14', '2022-11-14', '2022-12-21',...
    '2022-12-22', '2023-01-03', '2023-01-03', '2023-01-23', '2023-01-25'});
secondGroup = {'M045', 'M046', 'M047', 'M048'};

groupOne = ~ismember(ratListAll, secondGroup);
groupTwo = ismember(ratListAll, secondGroup);

[S, sDates,...
    sensitivityProcessed,...
    sensitivityAligned] = deal(cell(size(ratListAll)));

for rr = 1:length(ratListAll)
    ratname = ratListAll{rr};

    s = load([datadir 'S_Structs' filesep 'SStruct_' ratname]);
    S{rr} = s.S;

    sDates{rr} = cellfun(@(s) datetime(s.SessionDate), S{rr}.pd);
    sensitivityRaw = nan(size(sDates{rr}));

    for dd = 1:length(sDates{rr})
        disp([rr dd])
        Sday = struct('pd', {S{rr}.pd(dd)},...
            'peh', {S{rr}.peh(dd)});
        data = parse_data_from_mysql(Sday, [], 'all');

        latencyByPrevVolPre =...
            latency_post20ul_by_previous_vol(data);
        sensitivityRaw(dd) =...
            latencyByPrevVolPre(1)-latencyByPrevVolPre(2);
    end

    sensitivityProcessed{rr} = smooth(sensitivityRaw, 10, 'rloess');

    i = find(sDates{rr}>=startDates(rr), 1);
    sensitivityAligned{rr} = sensitivityProcessed{rr}(i:end);
end

%%
[SPre, dataPre, SPost, dataPost] = deal(cell(size(ratListAll)));

for rr = 1:length(ratListAll)
    disp(rr)
    CPCueDates = sDates{rr}'>=startDates(rr);

    if ~all(CPCueDates)
        SPre{rr} =...
            struct('pd', {S{rr}.pd(~CPCueDates)},...
            'peh', {S{rr}.peh(~CPCueDates)});
        dataPre{rr} = parse_data_from_mysql(SPre{rr});
    end

    SPost{rr} =...
        struct('pd', {S{rr}.pd(CPCueDates)},...
        'peh', {S{rr}.peh(CPCueDates)}); %#ok<*AGROW>
    dataPost{rr} = parse_data_from_mysql(SPost{rr}, [], 'all');
end

%%

[rewITIPreMean, rewITIPreSEM,...
    rewITIPostMean, rewITIPostSEM] =...
    deal(nan(length(ratListAll), 5));

[latencyByPrevVolPreMean, latencyByPrevVolPreSEM,...
    latencyByPrevVolPostMean, latencyByPrevVolPostSEM] =...
    deal(nan(length(ratListAll), 2));

[betasPreMean, betasPreSEM,...
    betasPostMean, betasPostSEM,...
    wt_betasMean, wt_betasSEM] =...
    deal(nan(length(ratListAll), nback+1));

[rewITIPrePval, rewITIPostPval,...
    latencyByPreVolPrePval, latencyByPreVolPostPval] =...
    deal(nan(length(ratListAll), 1));

[ltomPre, htomPre, mtolPre, mtohPre,...
    ltomPost, htomPost, mtolPost, mtohPost] =...
    deal(nan(length(ratListAll), 2*twin+1));

for rr = 1:length(ratListAll)
    disp(rr)

    if ~ismember(ratListAll{rr}, secondGroup)
        % Trial initiation time by offered reward
        [rewITIPreMean(rr,:),...
            rewITIPreSEM(rr,:),...
            rewITIPrePval(rr)] =...
            iti_by_reward(dataPre{rr});
        % Trial initiation time on 20 uL trials by previous offer
        [latencyByPrevVolPreMean(rr,:),...
            latencyByPrevVolPreSEM(rr,:),...
            latencyByPreVolPrePval(rr)] =...
            latency_post20ul_by_previous_vol(dataPre{rr});

        % Regression
        [betasPreMean(rr,:), ~, betasPreSEM(rr,:)] =...
            regress_latency_vs_rew(dataPre{rr}, nback, false, true);
        [wt_betasMean(rr,:), ~, wt_betasSEM(rr,:)] =...
            regress_wt_vs_rew(dataPre{rr}, nback, false, true);

        % Block Transition
        [ltomPre(rr,:), htomPre(rr,:), mtolPre(rr,:), mtohPre(rr,:)] =...
            block_dynamics_latency(dataPre{rr}, twin, smoothfactor);
    end

    [rewITIPostMean(rr,:),...
        rewITIPostSEM(rr,:),...
        rewITIPostPval(rr)] =...
        iti_by_reward(dataPost{rr});

    [latencyByPrevVolPostMean(rr,:),...
        latencyByPrevVolPostSEM(rr,:),...
        latencyByPreVolPostPval(rr)] =...
        latency_post20ul_by_previous_vol(dataPost{rr});

    [betasPostMean(rr,:), ~, betasPostSEM(rr,:)] =...
        regress_latency_vs_rew(dataPost{rr}, nback, false, true);

    [ltomPost(rr,:), htomPost(rr,:), mtolPost(rr,:), mtohPost(rr,:)] =...
        block_dynamics_latency(dataPost{rr}, twin, smoothfactor);
end

%%
clear l

npre = 12;
npost = 16;

figure
%--------------------------------------------------------------------------
% ITI by offered reward, pre-trained
%--------------------------------------------------------------------------
subplot(2, 2, 1); hold on
shadedErrorBar(1:5,...
    mean(rewITIPreMean(groupOne,:)),...
    std(rewITIPreMean(groupOne,:))./sqrt(npre))
shadedErrorBar(1:5,...
    mean(rewITIPostMean(groupOne,:), 'omitnan'),...
    std(rewITIPostMean(groupOne,:))./sqrt(npre),...
    lineprops='-m')

xlim([0.5 5.5])
ylim([-0.17 0.14])

title('pre-trained rats')

xticks(1:5)
xticklabels({'5', '10', '20', '40', '80'})
ylabel('Trial initiation time (z-score)')
xlabel('Offered volume')

%--------------------------------------------------------------------------
% ITI by reward, from scratch
%--------------------------------------------------------------------------
subplot(2, 2, 2); hold on
shadedErrorBar(1:5,...
    mean(rewITIPostMean(groupTwo,:), 'omitnan'),...
    std(rewITIPostMean(groupTwo,:))./sqrt(sum(groupTwo)),...
    lineprops='-m')
xlim([0.5 5.5])
ylim([-0.17 0.14])

title('only CPCue rats')

xticks(1:5)
xticklabels({'5', '10', '20', '40', '80'})

xlabel('Offered volume')

%--------------------------------------------------------------------------
% Conditional trial initiation time
%--------------------------------------------------------------------------
subplot(2, 2, 3); hold on
shadedErrorBar(1:2,...
    mean(latencyByPrevVolPreMean, 'omitnan'),...
    std(latencyByPrevVolPreMean, 'omitnan')./sqrt(npre),...
    lineprops='-k')
shadedErrorBar(1:2,...
    mean(latencyByPrevVolPostMean, 'omitnan'),...
    std(latencyByPrevVolPostMean, 'omitnan')./sqrt(npre),...
    lineprops='-m')

xticks(1:2)
xticklabels({'<20', '>20'})
xlim([0.5 2.5])

%--------------------------------------------------------------------------
% Block transitions pre- and post-CPCue
%--------------------------------------------------------------------------
subplot(2, 2, 4); hold on
shadedErrorBar(-twin:twin,...
    mean(ltomPre, 'omitnan'),...
    std(ltomPre, 'omitnan')./sqrt(npre),...
    lineprops='-b')
shadedErrorBar(-twin:twin,...
    mean(ltomPost, 'omitnan'),...
    std(ltomPost, 'omitnan')./sqrt(npost),...
    lineprops='--b')

shadedErrorBar(-twin:twin,...
    mean(htomPre, 'omitnan'),...
    std(htomPre, 'omitnan')./sqrt(npre),...
    lineprops='-r')
shadedErrorBar(-twin:twin,...
    mean(htomPost, 'omitnan'),...
    std(htomPost, 'omitnan')./sqrt(npost),...
    lineprops='--r')

set(gcf, Position=[1197 553 484 392])
