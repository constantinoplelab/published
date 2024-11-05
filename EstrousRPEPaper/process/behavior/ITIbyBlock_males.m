function ITIbyBlock = ITIbyBlock_males(m_ratlist, RatBehaviorData)

%Create variables
nback = 7;
prctile_cutoff = 98;
mrats = length(m_ratlist);

%create functions
block_avg = @(x, blocks) arrayfun(@(y) mean(x(blocks==y), 'omitnan'), 1:3);
block_err = @(x, blocks) arrayfun(@(y)...
    std(x(blocks==y), 'omitnan')./sqrt(sum(~isnan(x(blocks==y)))), 1:3);
rewards = [5 10 20 40 80];
reward_avg = @(x, rew) arrayfun(@(y) mean(x(rew==y), 'omitnan'), rewards);

%Get data from each rat
ITIbyBlock = [];
for rat = 1:mrats

    ratname = m_ratlist{rat};
    disp([ratname ' ' num2str(rat) ' out of ' num2str(mrats)])

    % Load data and process
    A = RatBehaviorData.(ratname);
    [~, A.reward] = convertreward(A.reward); %using standard reward amounts
    A.ITI(cumsum([0; A.ntrials(1:end-1)])+1) = NaN;
    A.ITI(A.ITI>prctile(A.ITI, prctile_cutoff)) = NaN;
    l = A.ITI;

    %%Calculate regression of latency as a function of trial
    %%number (in the future could split by thirst level?)
    %initialize matrix for latencies
    numsess = length(A.ntrials);
    sess = nan(numsess, max(A.ntrials));
    ctr = 1;
    xvec = (1:max(A.ntrials))';
    s = []; x = [];
    for jk = 1:numsess
        sess(jk,1:A.ntrials(jk)) = A.ITI(ctr:ctr+A.ntrials(jk)-1);
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

    %Detrend ITIs by session, accounting for satiety over the session and
    %remove outliers
    lat_cell = cell(numsess,1);
    lat_detrended_cell = cell(numsess,1);
    blocks_over_sessions = cell(numsess, 1);
    ctr = 1;
    for j = 1:numsess
        %Remove average trend over sessions
        numtrials = A.ntrials(j);
        ITIs_processed = A.ITI(ctr:ctr+numtrials-1);
        lat_cell{j} = ITIs_processed;
        lat_detrended_cell{j} = ITIs_processed - regline(1:numtrials)';

        A.ITI(ctr:ctr+numtrials-1) = lat_detrended_cell{j};
        blocks_over_sessions{j} = A.block(ctr:ctr+numtrials-1);
        ctr = ctr+numtrials;
    end
    latmat = cell2mat(lat_cell);
    lat_detrended_mat = cell2mat(lat_detrended_cell);
    blocktypes = cell2mat(blocks_over_sessions);

    %Get trial initiation times for low and high blocks
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
    [ratbetas, rattaus, ~, ratparamAs, ~, bestfit] =...
        regress_latency_vs_rew(B, nback, 0, 0, 0, 0); %use all mixed block trials
    [ratbetas_detrended, rattaus_detrended, ~, ratparamAs_detrended] =...
        regress_latency_vs_rew(B, nback, 0, 0, 0, 0);...% use all mixed block trials

    %get detrended ITI by previous reward for example plot inset
    mixedblock_rewards = B.reward(B.block==1);
    ITI_by_prevrew = reward_avg(B.ITI(B.block==1), [0; mixedblock_rewards(1:end-1)]);

    ITIbyBlock.('Male').high_ITI{rat} = avgITIbyblock(1, 2);
    ITIbyBlock.('Male').low_ITI{rat} = avgITIbyblock(1, 3);
    ITIbyBlock.('Male').high_ITI_err{rat} = errITIbyblock(1, 2);
    ITIbyBlock.('Male').low_ITI_err{rat} = errITIbyblock(1, 3);
    ITIbyBlock.('Male').high_ITI_det{rat} = avgITIbyblock_detrended(1, 2);
    ITIbyBlock.('Male').low_ITI_det{rat} = avgITIbyblock_detrended(1, 3);
    ITIbyBlock.('Male').high_ITI_det_err{rat} = errITIbyblock_detrended(1, 2);
    ITIbyBlock.('Male').low_ITI_det_err{rat} = errITIbyblock_detrended(1, 3);
    ITIbyBlock.('Male').ranksum_pval_det(rat) = ranksum_pval_det;
    ITIbyBlock.('Male').deltas(rat) = delta;
    ITIbyBlock.('Male').deltas_det(rat) = delta_detrended;
    ITIbyBlock.('Male').betas(rat,:) = ratbetas';
    ITIbyBlock.('Male').taus(rat,1) = rattaus;
    ITIbyBlock.('Male').As(rat,1) = ratparamAs;
    ITIbyBlock.('Male').ITI_by_prevrew(rat,:) = ITI_by_prevrew;
    ITIbyBlock.('Male').bestfit(rat,:) = bestfit';
    ITIbyBlock.('Male').betas_detrended(rat,:) = ratbetas_detrended';
    ITIbyBlock.('Male').taus_detrended(rat,1) = rattaus_detrended;
    ITIbyBlock.('Male').paramAs_detrended(rat,1) = ratparamAs_detrended;
end

end
