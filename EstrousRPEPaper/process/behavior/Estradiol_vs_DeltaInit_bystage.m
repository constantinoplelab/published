function [delta_rats, E2_rats, ratlist_elisa] =...
    Estradiol_vs_DeltaInit_bystage(ratlist,...
    SerumTable, SerumBehData, detrend_arg)
%relate change in estradiol expression to change in delta initiation times, ie "behavioral sensitivity" (pro -
%di) only for sessions with blood samples
%calculates mean concentration for repeated dates

%only keep full samples (50 ul)
fullsample = SerumTable.Volume==50;
SerumTable = SerumTable(fullsample, :);

%calculate mean high and low detrended initiation times for each session
prctile_cutoff = 99;
ratlist_elisa = unique(SerumTable.Rat);
ratlist_elisa = ratlist_elisa(ismember(ratlist_elisa, ratlist)); %only keep rats for this experiment
low_rats = NaN(length(ratlist_elisa), 2);
high_rats = NaN(length(ratlist_elisa), 2);
E2_rats = NaN(length(ratlist_elisa), 2);
delta_rats = NaN(length(ratlist_elisa), 2);
cycle = {'Proestrus', 'Diestrus'};
for rat = 1:length(ratlist_elisa)

    if contains(ratlist_elisa{rat}, 'June')
        fullname = ratlist_elisa{rat};
        disp(fullname)
    else
        fullname = ratlist_elisa{rat};
        ratname = ratlist_elisa{rat};
        disp(ratname)
    end

    %load rat's behavioral data
    stagedS = SerumBehData.(ratname);

    %remove empty stage fields
    emptystageidx = cellfun(@(x) ~isfield(x, 'Stage'), stagedS.pd);
    stagedS.pd(emptystageidx) = [];
    stagedS.peh(emptystageidx) = [];

    for e = 1:length(cycle)

        %subset to this stage
        stageidx = cellfun(@(x) logical(sum(strcmp(x.Stage,...
            cycle{e}))), stagedS.pd); %find where it matches cycle stage(s)
        stageS.pd = stagedS.pd(stageidx);
        stageS.peh = stagedS.peh(stageidx);

        %subset to sampled sessions only
        serum_sessiondates =...
            SerumTable.SessionDate(strcmp(SerumTable.Rat, fullname)...
            & strcmp(SerumTable.Stage, cycle{e}));
        serum_sessiondates_unique = unique(serum_sessiondates);
        behsess_withserum_idx = find(cellfun(@(x)...
            logical(ismember(x.SessionDate, serum_sessiondates_unique)),...
            stageS.pd));
        disp([num2str(length(serum_sessiondates_unique))...
            ' sessions captured'])
        sampledS = stageS;
        sampledS.pd = sampledS.pd(behsess_withserum_idx);
        sampledS.peh = sampledS.peh(behsess_withserum_idx);

        %convert to A struct to be able to detrend
        A = parse_data_from_mysql(sampledS, 1, 1);
        [~, A.reward] = convertreward(A.reward); %using standard reward amounts
        A.ITI(cumsum([0; A.ntrials(1:end-1)])+1) = NaN;
        A.ITI(A.ITI>prctile(A.ITI, prctile_cutoff)) = NaN;
        l = A.ITI;

        %%we want to regress latency against trial number, so we can subtract
        %%it out later.
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
        bad = find(isnan(s));
        s(bad) = [];
        x(bad) = [];
        X = [x, ones(length(x),1)];
        [b] = regress(s,X);
        regline = (b(1).*xvec + b(2))'; %average latency by trial number

        %calculate ITI vector across sessions and NaN first trial
        numsess = length(sampledS.pd);
        [ITIs_over_sessions, ~, ~] = getITI(sampledS); %includes post-violations
        for j = 1:numsess
            ITIs_thissess = ITIs_over_sessions{j};
            ITIs_NaNfirsttrial = [NaN; ITIs_thissess(2:end)];
            ITIs_over_sessions{j} = ITIs_NaNfirsttrial;
        end
        ITIs = cell2mat(ITIs_over_sessions');

        %Find ITI threshold to cutoff
        lat_thresh = prctile(ITIs,prctile_cutoff);

        %calculate for each session
        low_sess = cell(length(serum_sessiondates_unique), 1);
        high_sess = cell(length(serum_sessiondates_unique), 1);
        for sess = 1:length(serum_sessiondates_unique)

            %detrend across session
            ITIs_processed = ITIs_over_sessions{sess};
            ITIs_processed(ITIs_processed>lat_thresh) = NaN;
            if detrend_arg == 1
                lat = ITIs_processed - regline(1:length(ITIs_processed))';
            else
                lat = ITIs_processed;
            end

            blocks = sampledS.pd{sess}.Block;

            %calculate detrended mean low and hi and plot
            low_sess{sess} = lat(blocks==3);
            high_sess{sess} = lat(blocks==2);

        end
        low_mat = cell2mat(low_sess);
        high_mat = cell2mat(high_sess);

        low_rats(rat, e) = mean(low_mat, 'omitnan');
        high_rats(rat, e) = mean(high_mat, 'omitnan');
        delta_rats(rat, e) = low_rats(rat, e) - high_rats(rat, e);
        E2_rats(rat, e) =...
            mean(SerumTable.EstradiolConc(strcmp(SerumTable.Rat, fullname)...
            & strcmp(SerumTable.Stage, cycle{e})), 'omitnan');

    end
end

end