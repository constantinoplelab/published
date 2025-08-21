function ITIbyBlock = ITIbyBlock_females(f_ratlist, RatBehaviorData)

%Create variables
nback = 7;
prctile_cutoff = 98;
frats = length(f_ratlist);

%Create functions
block_avg = @(x, blocks) arrayfun(@(y) mean(x(blocks==y), 'omitnan'), 1:3);
block_err = @(x, blocks) arrayfun(@(y)...
    std(x(blocks==y), 'omitnan')./sqrt(sum(~isnan(x(blocks==y)))), 1:3);

%Get data from each rat
ITIbyBlock = [];
for rat = 1:frats
    ratname = f_ratlist{rat};
    disp([ratname ' ' num2str(rat) ' out of ' num2str(frats)])

    S = RatBehaviorData.(ratname);
    S = S.S;

    %remove empty stage fields
    emptystageidx = cellfun(@(x) ~isfield(x, 'Stage'), S.pd);
    stagedS = S;
    stagedS.pd(emptystageidx) = [];
    stagedS.peh(emptystageidx) = [];

    numsess = length(stagedS.pd);

    %%Calculate regression of latency as a function of trial
    %%number for sessions from this stage
    A = parse_data_from_mysql(stagedS, 1, 1);
    [~, A.reward] = convertreward(A.reward); %using standard reward amounts
    A.ITI(cumsum([0; A.ntrials(1:end-1)])+1) = NaN;
    A.ITI(A.ITI>prctile(A.ITI, prctile_cutoff)) = NaN;
    l = A.ITI;

    %initialize matrix for latencies
    sess = nan(length(A.ntrials), max(A.ntrials));
    ctr = 1;
    xvec = (1:max(A.ntrials))';
    s = []; x = [];
    for jk = 1:length(A.ntrials)
        sess(jk,1:A.ntrials(jk)) = l(ctr:ctr+A.ntrials(jk)-1);
        x = [x; xvec];
        s = [s; sess(jk,:)'];
        ctr = ctr+A.ntrials(jk);
    end
    %%we want to regress latency against trial number, so we can subtract
    %%it out later.
    bad = find(isnan(s));
    s(bad) = [];
    x(bad) = [];
    X = [x, ones(length(x),1)];
    [b] = regress(s,X);
    regline = (b(1).*xvec + b(2))'; %average latency by trial number

    %calculate ITI vector across sessions and NaN first trial
    [ITIs_over_sessions, ~, ~] = getITI(stagedS); %includes post-violations
    for j = 1:numsess
        ITIs_thissess = ITIs_over_sessions{j};
        ITIs_NaNfirsttrial = [NaN; ITIs_thissess(2:end)];
        ITIs_over_sessions{j} = ITIs_NaNfirsttrial;
    end
    ITIs = cell2mat(ITIs_over_sessions');

    %Find ITI threshold to cutoff
    lat_thresh = prctile(ITIs,prctile_cutoff);

    %Detrend ITIs, accounting for satiety over the session and
    %remove outliers
    lat_cell = cell(numsess,1);
    lat_detrended_cell = cell(numsess,1);
    blocks_over_sessions = cell(numsess, 1);
    for j = 1:numsess
        %remove outliers and effect of satiety over session
        ITIs_processed = ITIs_over_sessions{j};
        ITIs_processed(ITIs_processed>lat_thresh) = NaN;
        lat_cell{j} = ITIs_processed;
        lat_detrended_cell{j} = ITIs_processed - regline(1:length(ITIs_processed))';
        blocks_over_sessions{j} = stagedS.pd{j}.Block;
    end
    latmat = cell2mat(lat_cell);
    lat_detrended_mat = cell2mat(lat_detrended_cell);
    blocktypes = cell2mat(blocks_over_sessions);

    %Get trial initiation times for mixed, low and high blocks
    avgITIbyblock = block_avg(latmat, blocktypes);
    errITIbyblock = block_err(latmat, blocktypes);
    delta = mean(latmat(blocktypes==3), 'omitnan') -...
        mean(latmat(blocktypes==2), 'omitnan');

    %Get trial initiation times for low and high blocks
    %(detrended)
    avgITIbyblock_detrended = block_avg(lat_detrended_mat, blocktypes);
    errITIbyblock_detrended = block_err(lat_detrended_mat, blocktypes);
    delta_detrended = mean(lat_detrended_mat(blocktypes==3), 'omitnan') -...
        mean(lat_detrended_mat(blocktypes==2), 'omitnan');

    %Determine whether low and high block initiation times are
    %significantly different
    if sum(blocktypes==2) > 0 && sum(blocktypes==3) > 0
        ranksum_pval_det = ranksum(lat_detrended_mat(blocktypes==2),...
            lat_detrended_mat(blocktypes==3));
    else
        ranksum_pval_det = NaN;
    end
    
    %Subtract latency as a function of trial number
    ctr = 1;
    B = A;
    for jk = 1:length(B.ntrials)
        %Remove average trend over sessions
        B.ITI(ctr:ctr+B.ntrials(jk)-1) =...
            l(ctr:ctr+B.ntrials(jk)-1) - regline(1:B.ntrials(jk))';
        ctr = ctr+B.ntrials(jk);
    end

    %Regress (no detrending)
    [ratbetas, rattaus, ~, ratparamAs] =...
        regress_latency_vs_rew(A, nback, 0, 0, 0, 0); %, usethese
    [ratbetas_detrended, rattaus_detrended, ~, ratparamAs_detrended] =...
        regress_latency_vs_rew(B, nback, 0, 0, 0, 0);...% use all mixed block trials
        
    ITIbyBlock.('Female').mixed_ITI{rat} = avgITIbyblock(1, 1);
    ITIbyBlock.('Female').high_ITI{rat} = avgITIbyblock(1, 2);
    ITIbyBlock.('Female').low_ITI{rat} = avgITIbyblock(1, 3);
    ITIbyBlock.('Female').mixed_ITI_err{rat} = errITIbyblock(1, 1);
    ITIbyBlock.('Female').high_ITI_err{rat} = errITIbyblock(1, 2);
    ITIbyBlock.('Female').low_ITI_err{rat} = errITIbyblock(1, 3);
    ITIbyBlock.('Female').mixed_ITI_det{rat} = avgITIbyblock_detrended(1, 1);
    ITIbyBlock.('Female').high_ITI_det{rat} = avgITIbyblock_detrended(1, 2);
    ITIbyBlock.('Female').low_ITI_det{rat} = avgITIbyblock_detrended(1, 3);
    ITIbyBlock.('Female').mixed_ITI_det_err{rat} = errITIbyblock_detrended(1, 1);
    ITIbyBlock.('Female').high_ITI_det_err{rat} = errITIbyblock_detrended(1, 2);
    ITIbyBlock.('Female').low_ITI_det_err{rat} = errITIbyblock_detrended(1, 3);
    ITIbyBlock.('Female').ranksum_pval_det(rat) = ranksum_pval_det;
    ITIbyBlock.('Female').deltas(rat) = delta;
    ITIbyBlock.('Female').deltas_det(rat) = delta_detrended;
    ITIbyBlock.('Female').betas(rat,:) = ratbetas';
    ITIbyBlock.('Female').taus(rat,1) = rattaus;
    ITIbyBlock.('Female').As(rat,1) = ratparamAs;
    ITIbyBlock.('Female').betas_detrended(rat,:) = ratbetas_detrended';
    ITIbyBlock.('Female').taus_detrended(rat,1) = rattaus_detrended;
    ITIbyBlock.('Female').paramAs_detrended(rat,1) = ratparamAs_detrended;

end
