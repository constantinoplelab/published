function ITIbyBlock = ITIbyBlock_dox_stages(ratlist,...
    lentirats, DoxBeh)
%input:
%ratlist = dox rats (e.g. {'G067', 'G069', 'G076'})
%project = shEsr1 or control
%% GET DATA

%settings
% twin = 18;
twin = 40;
smoothfactor = 15;
nback = 7;
prctile_cutoff = 99;
cycle = {'Proestrus', 'Diestrus'};
cyclenames = cycle;

%create functions
block_avg = @(x, blocks) arrayfun(@(y) mean(x(blocks==y), 'omitnan'), 1:3);
block_err = @(x, blocks) arrayfun(@(y)...
    std(x(blocks==y), 'omitnan')./sqrt(sum(~isnan(x(blocks==y)))), 1:3);
numsess_rats = cell(length(ratlist), 1);

for rat = 1:length(ratlist)

    ratname = char(ratlist(rat));
    disp(ratname)

    if ismember(ratname, lentirats) %wait three weeks after lentiviral expression
        dayspostdox = 21;
    else
        dayspostdox = 1;
    end

    S = DoxBeh.(ratname);

    %assign dox states to rats
    if strcmp(ratname, {'G079'})
        dox_states = {'predox', 'duringdox', 'postdox', 'sugar', 'doxnosugar'};
    elseif strcmp(ratname, {'G086'}) %didn't cycle while on doxnosugar
        dox_states = {'predox', 'duringdox', 'doxnosugar'};
    elseif ismember(ratname, {'G083','G100','G101','G110','G133', 'G134', 'G135', 'G136', 'G141'}) %never taken off of dox, G083 had sugar in dox, G133 got lenti-shEsr1
        dox_states = {'predox', 'duringdox'};
    elseif ismember(ratname, {'G091','G093','G094','G096','G097','G098','G113','G115'})
        dox_states = {'predox', 'duringdox', 'end_date'};
    elseif ismember(ratname, {'G108','G109','G120'})
        dox_states = {'predox', 'dox_controlrat_start', 'dox_controlrat_end'};
    elseif ismember(ratname, {'G112'})
        dox_states = {'predox', 'dox_controlrat_start', 'dox_controlrat_end', 'dox_controlrat_startagain'};
    elseif ismember(ratname, {'G092','G099','G111','G114','G116'}) % 99 did control sessions (ie no shRNA) on dox and then got shEsr2 injection, do not do 'G101' b/c only got control chow before surgery
        dox_states = {'predox', 'dox_controlrat_start', 'dox_controlrat_end'};
    end

    doxstatenames = {'predox','duringdox'};
    length_dox_states = 2;
    datesname = cell(length_dox_states,1);
    numsess_by_doxstate = zeros(length(cycle), length_dox_states);

    for d = 1:length_dox_states

        doxS = S.doxS;
        doxdates = get_dox_dates(ratlist(rat));
        disp(doxstatenames{d})

        if strcmp(doxstatenames(d), 'predox') %DOX STATE 1
            if  ~ismember('dox_controlrat_start', dox_states) %dox rats
                dateindex = cellfun(@(x) datetime(x.SessionDate)...
                    < datetime(doxdates.start_date), doxS.pd);
                datesname{d} = ['before dox ' doxdates.start_date{1}];
            elseif ismember('dox_controlrat_start', dox_states) %control rats
                dateindex = cellfun(@(x) datetime(x.SessionDate)...
                    < datetime(doxdates.dox_controlrat_start), doxS.pd);
                datesname{d} = ['before dox ' doxdates.dox_controlrat_start{1}];
            end
        elseif strcmp(doxstatenames(d), 'duringdox') %DOX STATE 2
            if ismember(ratname, lentirats)
                effectivestart = datetime(doxdates.start_date)+dayspostdox;
                dateindex = cellfun(@(x) datetime(x.SessionDate)...
                    > effectivestart, doxS.pd);
            else
                dox_effective_startdate = datetime(doxdates.start_date)+dayspostdox;
                if ~isempty(doxdates.end_date)
                    dateindex = cellfun(@(x) datetime(x.SessionDate)...
                        > dox_effective_startdate & datetime(x.SessionDate)...
                        < datetime(doxdates.end_date), doxS.pd);
                elseif isempty(doxdates.end_date)
                    dateindex = cellfun(@(x) datetime(x.SessionDate)...
                        > dox_effective_startdate, doxS.pd);
                end
                datesname{d} = ['after ' num2str(dayspostdox) ' days into dox (' ...
                    datestr(dox_effective_startdate) ') and before being taken off('...
                    datestr(doxdates.end_date) ')'];
            end
        end

        doxS.pd = doxS.pd(dateindex);
        doxS.peh = doxS.peh(dateindex);
        %remove empty stage fields
        emptystageidx = cellfun(@(x) ~isfield(x, 'Stage'), doxS.pd);
        doxS.pd(emptystageidx) = [];
        doxS.peh(emptystageidx) = [];

        numsess = length(doxS.pd);

        if numsess > 0

            %%Calculate regression of latency as a function of trial
            %%number for sessions from this stage
            SforA.pd = doxS.pd;
            SforA.peh = doxS.peh;
            A = parse_data_from_mysql(SforA, 1, 1);
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

            %get ITIs per session
            [ITIs_over_sessions, ~, ~] = getITI(doxS);
            %NaN first trial
            for j = 1:numsess
                ITIs_thissess = ITIs_over_sessions{j};
                ITIs_NaNfirsttrial = [NaN; ITIs_thissess(2:end)];
                ITIs_over_sessions{j} = ITIs_NaNfirsttrial;
            end
            ITIs = cell2mat(ITIs_over_sessions');

            %Find ITI threshold to cutoff
            lat_thresh = prctile(ITIs,prctile_cutoff);

            %Detrend ITIs, accounting for satiety over the session and
            %remove outliers, get thirst and motivation metrics
            lat_cell = cell(numsess,1);
            lat_detrended_cell = cell(numsess,1);
            blocks_over_sessions = cell(numsess, 1);
            numtrialscorrected = NaN(numsess, 1);
            volumecorrected = NaN(numsess, 1);
            for j = 1:numsess
                %remove outliers and effect of satiety over session
                ITIs_processed = ITIs_over_sessions{j};
                ITIs_processed(ITIs_processed>lat_thresh) = NaN;
                lat_cell{j} = ITIs_processed;
                lat_detrended_cell{j} = ITIs_processed - regline(1:length(ITIs_processed))';
                blocks_over_sessions{j} = doxS.pd{j}.Block;

                %get thirst and motivation metrics
                %avg number of trials per session
                numtrials = length(doxS.pd{j}.hits);
                %volume consumed
                vol = sum(doxS.pd{j}.RewardAmount(doxS.pd{j}.hits));
                RewardRateSess = DoxBeh.(ratname).RewardRateStruct.pd{j};
                if (RewardRateSess.TrialEnd(end) -...
                        RewardRateSess.TrialStart(1)) < 10800 %exclude sessions longer than 3 hours
                    duration = RewardRateSess.TrialEnd(end) -...
                        RewardRateSess.TrialStart(1);
                    numtrialscorrected(j) = numtrials./duration;
                    volumecorrected(j) = vol./duration;
                else
                    numtrialscorrected(j) = NaN;
                    volumecorrected(j) = NaN;
                end
            end
            lat_detrended_mat = cell2mat(lat_detrended_cell);
            blocktypes = cell2mat(blocks_over_sessions);

            ITIbyBlock.(doxstatenames{d}).NumTrials(rat) = mean(numtrials, 'omitnan');
            ITIbyBlock.(doxstatenames{d}).TrialsOverDuration(rat) = mean(numtrialscorrected, 'omitnan'); %divide volume and number of trials by session length (separately)
            ITIbyBlock.(doxstatenames{d}).VolumeOverDuration(rat) = mean(volumecorrected, 'omitnan'); %divide volume and number of trials by session length (separately)

            %Get trial initiation times for low and high blocks
            %(detrended)
            avgITIbyblock_detrended = block_avg(lat_detrended_mat, blocktypes);
            errITIbyblock_detrended = block_err(lat_detrended_mat, blocktypes);
            delta_detrended = mean(lat_detrended_mat(blocktypes==3), 'omitnan') -...
                mean(lat_detrended_mat(blocktypes==2), 'omitnan');

            ITIbyBlock.(doxstatenames{d}).high(rat) = avgITIbyblock_detrended(1, 2);
            ITIbyBlock.(doxstatenames{d}).low(rat) = avgITIbyblock_detrended(1, 3);
            ITIbyBlock.(doxstatenames{d}).high_err(rat) = errITIbyblock_detrended(1, 2);
            ITIbyBlock.(doxstatenames{d}).low_err(rat) = errITIbyblock_detrended(1, 3);
            ITIbyBlock.(doxstatenames{d}).delta(rat) = delta_detrended;

            %Regress latency vs. reward
            % First, subtract latency as a function of trial number
            ctr = 1;
            B = A;
            for jk = 1:length(B.ntrials)
                %Remove average trend over sessions
                B.ITI(ctr:ctr+B.ntrials(jk)-1) =...
                    l(ctr:ctr+B.ntrials(jk)-1) - regline(1:B.ntrials(jk))';
                ctr = ctr+B.ntrials(jk);
            end
            [ITIbyBlock.(doxstatenames{d}).betas(rat, :),~, ~] =...
                regress_latency_vs_rew(B, nback, 0, 0, 0, 0); %detrended

            %look at latency dynamics
            [ltom, htom, mtol, mtoh] = ... %decide in here whether or not to z-score, causal
                block_dynamics_latency(A, twin, smoothfactor, 0, 1, 0, 1);
            ITIbyBlock.(doxstatenames{d}).ltom(rat,:) = ltom;
            ITIbyBlock.(doxstatenames{d}).htom(rat,:) = htom;
            ITIbyBlock.(doxstatenames{d}).mtol(rat,:) = mtol;
            ITIbyBlock.(doxstatenames{d}).mtoh(rat,:) = mtoh;

            ITIbyBlock.numsess_doxstate(rat, d) = numsess;

            %% Separate by stage group
            if ~isempty(doxS.pd)

                %remove empty stage fields, in case you didn't before
                emptystageidx = cellfun(@(x) ~isfield(x, 'Stage'), doxS.pd);
                doxS.pd(emptystageidx) = [];
                doxS.peh(emptystageidx) = [];
                RewardRateStruct = DoxBeh.(ratname).RewardRateStruct;
                RewardRateStruct.pd(emptystageidx) = [];

                for e = 1:length(cycle)

                    stageidx = cellfun(@(x) logical(sum(strcmp(x.Stage,...
                        cycle{e}))), doxS.pd); %find where it matches cycle stage(s)
                    stageS.pd = doxS.pd(stageidx);
                    stageS.peh = doxS.peh(stageidx);
                    stageRewardRateStruct.pd = RewardRateStruct.pd(stageidx);
                    disp([cycle{e} ': ' num2str(sum(stageidx)) ' sessions'])

                    if ~isempty(stageS.pd)

                        numsess_by_doxstate(e, d) = length(stageS.pd);

                        %%Calculate regression of latency as a function of trial
                        %%number for sessions from this stage
                        A = parse_data_from_mysql(stageS, 1, 1);
                        [~, A.reward] = convertreward(A.reward); %using standard
                        B = A;
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

                        %calculate ITI vector across sessions
                        [ITIs_over_sessions, ~, ~] = getITI(stageS);

                        %NaN first trial
                        for j = 1:numsess_by_doxstate(e, d)
                            ITIs_thissess = ITIs_over_sessions{j};
                            ITIs_NaNfirsttrial = [NaN; ITIs_thissess(2:end)];
                            ITIs_over_sessions{j} = ITIs_NaNfirsttrial;
                        end
                        ITIs_processed = cell2mat(ITIs_over_sessions');
                        lat_thresh = prctile(ITIs_processed,prctile_cutoff);
                        % figure; histogram(ITIs_processed); xline(lat_thresh, '--r')

                        %Detrend ITIs, accounting for satiety over the session
                        numsess = length(stageS.pd);
                        lat_cell = cell(numsess,1);
                        blocks_over_sessions = cell(numsess, 1);
                        high_sess_avg = NaN(numsess, 1);
                        low_sess_avg = NaN(numsess, 1);
                        numtrialscorrected =  NaN(numsess, 1);
                        volumecorrected =  NaN(numsess, 1);
                        for j = 1:numsess
                            %remove outliers and effect of satiety over session
                            ITIs_processed = ITIs_over_sessions{j};
                            ITIs_processed(ITIs_processed>lat_thresh) = NaN;
                            lat_cell{j} = ITIs_processed - regline(1:length(ITIs_processed))';
                            these_lat = ITIs_processed - regline(1:length(ITIs_processed))';
                            blocks_over_sessions{j} = stageS.pd{j}.Block;
                            these_blk = stageS.pd{j}.Block;

                            %get session averages
                            high_sess_avg(j, 1) = mean(these_lat(these_blk==2), 'omitnan');
                            low_sess_avg(j, 1) = mean(these_lat(these_blk==3), 'omitnan');

                            %get thirst and motivation metrics
                            %avg number of trials per session
                            numtrials = length(stageS.pd{j}.hits);
                            %volume consumed
                            vol = sum(stageS.pd{j}.RewardAmount(stageS.pd{j}.hits));
                            if (stageRewardRateStruct.pd{j}.TrialEnd(end) -...
                                    stageRewardRateStruct.pd{j}.TrialStart(1)) < 10800 %exclude sessions longer than 3 hours
                                duration = stageRewardRateStruct.pd{j}.TrialEnd(end) -...
                                    stageRewardRateStruct.pd{j}.TrialStart(1);
                                numtrialscorrected(j) = numtrials./duration;
                                volumecorrected(j) = vol./duration;
                            else
                                numtrialscorrected(j) = NaN;
                                volumecorrected(j) = NaN;
                            end
                        end
                        latmat = cell2mat(lat_cell);
                        blocktypes = cell2mat(blocks_over_sessions);

                        %average by block
                        avgITIbyblock = block_avg(latmat, blocktypes);
                        %error by block
                        errITIbyblock = block_err(latmat, blocktypes);
                        %calculate delta
                        delta = mean(latmat(blocktypes==3), 'omitnan') -...
                            mean(latmat(blocktypes==2), 'omitnan');

                        %%Regress latency vs. reward
                        %Subtract latency as a function of trial number
                        ctr = 1;
                        for jk = 1:length(A.ntrials)
                            %Remove average trend over sessions
                            A.ITI(1:A.ntrials(jk), 1) =...
                                l(ctr:ctr+A.ntrials(jk)-1) - regline(1:A.ntrials(jk))';
                            ctr = ctr+A.ntrials(jk);
                        end
                        [ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).betas(rat, :),~, ~] =...
                            regress_latency_vs_rew(A, nback, 0, 0, 0);

                        %look at latency dynamics
                        [ltom, htom, mtol, mtoh] = ... %decide in here whether or not to z-score
                            block_dynamics_latency(B, twin, smoothfactor, 0, 1, 1); %nonpostvio

                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).ltom(rat,:) = ltom;
                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).htom(rat,:) = htom;
                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).mtol(rat,:) = mtol;
                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).mtoh(rat,:) = mtoh;
                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).delta(rat,:) = delta;
                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).high(rat) = avgITIbyblock(1, 2);
                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).high_sess{rat} = high_sess_avg;
                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).low(rat) = avgITIbyblock(1, 3);
                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).low_sess{rat} = low_sess_avg;
                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).high_err(rat) = errITIbyblock(1, 2);
                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).low_err(rat) = errITIbyblock(1, 3);
                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).TrialsOverDuration(rat) =...
                            mean(numtrialscorrected, 'omitnan'); %divide volume and number of trials by session length (separately)
                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).VolumeOverDuration(rat) =...
                            mean(volumecorrected, 'omitnan'); %divide volume and number of trials by session length (separately)

                    else
                        disp(['no ' cyclenames{e} ' sessions'])
                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).ltom(rat,:) = NaN(1, 2*twin+1);
                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).htom(rat,:) = NaN(1, 2*twin+1);
                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).mtol(rat,:) = NaN(1, 2*twin+1);
                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).mtoh(rat,:) = NaN(1, 2*twin+1);
                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).delta(rat,:) = NaN;
                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).betas(rat,:) = NaN(1, nback+1);
                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).high(rat) = NaN;
                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).low(rat) = NaN;
                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).high_err(rat) = NaN;
                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).low_err(rat) = NaN;
                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).VolumeOverDuration(rat) = NaN;
                        ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).TrialsOverDuration(rat) = NaN;
                    end
                end
            else
                for e = 1:length(cycle)
                    ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).ltom(rat,:) = NaN(1, 2*twin+1);
                    ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).htom(rat,:) = NaN(1, 2*twin+1);
                    ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).mtol(rat,:) = NaN(1, 2*twin+1);
                    ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).mtoh(rat,:) = NaN(1, 2*twin+1);
                    ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).delta(rat,:) = NaN;
                    ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).betas(rat,:) = NaN(1, nback+1);
                    ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).high(rat) = NaN;
                    ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).low(rat) = NaN;
                    ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).high_err(rat) = NaN;
                    ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).low_err(rat) = NaN;
                    ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).VolumeOverDuration(rat) = NaN;
                    ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).TrialsOverDuration(rat) = NaN;
                end
            end
            numsess_rats{rat, 1} = numsess_by_doxstate;
        else
            ITIbyBlock.(doxstatenames{d}).volumeconsumed(rat) = NaN;
            ITIbyBlock.(doxstatenames{d}).numtrials(rat) = NaN;
            numsess_rats{rat, 1} = NaN;
            for e = 1:length(cycle)
                ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).ltom(rat,:) = NaN(1, 2*twin+1);
                ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).htom(rat,:) = NaN(1, 2*twin+1);
                ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).mtol(rat,:) = NaN(1, 2*twin+1);
                ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).mtoh(rat,:) = NaN(1, 2*twin+1);
                ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).delta(rat,:) = NaN;
                ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).betas(rat,:) = NaN(1, nback+1);
                ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).high(rat) = NaN;
                ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).low(rat) = NaN;
                ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).high_err(rat) = NaN;
                ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).low_err(rat) = NaN;
                ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).VolumeOverDuration(rat) = NaN;
                ITIbyBlock.(doxstatenames{d}).(cyclenames{e}).TrialsOverDuration(rat) = NaN;
            end
            numsess_rats{rat, 1} = 0;
            ITIbyBlock.numsess_doxstate(rat, d) = 0;
            ITIbyBlock.(doxstatenames{d}).ltom(rat,:) = NaN(1, 2*twin+1);
            ITIbyBlock.(doxstatenames{d}).htom(rat,:) = NaN(1, 2*twin+1);
            ITIbyBlock.(doxstatenames{d}).mtol(rat,:) = NaN(1, 2*twin+1);
            ITIbyBlock.(doxstatenames{d}).mtoh(rat,:) = NaN(1, 2*twin+1);
            ITIbyBlock.(doxstatenames{d}).high(rat) = NaN;
            ITIbyBlock.(doxstatenames{d}).low(rat) = NaN;
            ITIbyBlock.(doxstatenames{d}).high_err(rat) = NaN;
            ITIbyBlock.(doxstatenames{d}).low_err(rat) = NaN;
            ITIbyBlock.(doxstatenames{d}).delta(rat) = NaN;
            ITIbyBlock.(doxstatenames{d}).betas(rat,:) = NaN(1, nback+1);
            ITIbyBlock.(doxstatenames{d}).VolumeOverDuration(rat) = NaN;
            ITIbyBlock.(doxstatenames{d}).TrialsOverDuration(rat) = NaN;
            ITIbyBlock.(doxstatenames{d}).NumTrials(rat) = NaN;
        end
    end

    ITIbyBlock.numsess_doxstate_stage{rat} = numsess_by_doxstate;

end

