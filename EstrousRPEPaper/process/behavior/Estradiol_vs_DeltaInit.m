function [x_good, y_good, stages_good,... %individual sessions
    xbin, ybin, ybiner] =... %binned data
    Estradiol_vs_DeltaInit(SerumTable, SerumBehData)
%correlate change in estradiol expression to changes in delta initiation times, ie "behavioral sensitivity"
% where each point is a session and rats are pooled
%calculates mean concentration for repeated dates

%only keep full samples (50 ul)
fullsample = SerumTable.Volume==50;
SerumTable = SerumTable(fullsample, :);

%calculate mean high and low detrended initiation times for each session
prctile_cutoff = 98;
disp(num2str(prctile_cutoff))
ratlist_elisa = unique(SerumTable.Rat);
delta_sess_overrats = [];
E2_sess_overrats = [];
stages_sess_overrats = [];
for rat = 1:length(ratlist_elisa)

    fullname = ratlist_elisa{rat};
    ratname = ratlist_elisa{rat};
    disp(ratname)

    %load rat's behavioral data
    stagedS = SerumBehData.(ratname);

    %remove empty stage fields
    emptystageidx = cellfun(@(x) ~isfield(x, 'Stage'), stagedS.pd);
    stagedS.pd(emptystageidx) = [];
    stagedS.peh(emptystageidx) = [];

    %subset to sampled sessions only
    serum_sessiondates = SerumTable.SessionDate(strcmp(SerumTable.Rat, fullname));
    serum_sessiondates_unique = unique(serum_sessiondates); %orders them as well
    behsess_withserum_idx = find(cellfun(@(x)...
        logical(ismember(x.SessionDate, serum_sessiondates_unique)),...
        stagedS.pd));
    disp([num2str(length(serum_sessiondates_unique)) ' sessions captured'])
    sampledS = stagedS;
    sampledS.pd = sampledS.pd(behsess_withserum_idx);
    sampledS.peh = sampledS.peh(behsess_withserum_idx);
    dates_beh = cellfun(@(x) x.SessionDate, sampledS.pd,...
        'UniformOutput', false);
    disp(dates_beh)

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
    disp(lat_thresh)

    %calculate delta for each session
    %detrend for each session
    delta_sess = cell(length(serum_sessiondates_unique), 1);
    stages_rat = cell(length(serum_sessiondates_unique), 1);
    for sess = 1:length(serum_sessiondates_unique)
            thissess.pd = sampledS.pd(sess);
            thissess.peh = sampledS.peh(sess);

            ITIs_processed = ITIs_over_sessions{sess};
            ITIs_processed(ITIs_processed>lat_thresh) = NaN; %remove outliers

            % %detrending
            % %detrend using the regression line fit to the mean of all
            % sessions.
            lat = ITIs_processed - regline(1:length(ITIs_processed))';

            blocks = thissess.pd{1}.Block;

            %calculate detrended delta
            delta_sess{sess} = (mean(lat(blocks==3), 'omitnan') -...
                mean(lat(blocks==2), 'omitnan'));

            stages_rat{sess} = thissess.pd{1}.Stage;
        
    end

    delta_sess_overrats = [delta_sess_overrats; cell2mat(delta_sess)];
    ratTable = SerumTable(strcmp(SerumTable.Rat, fullname), :); 
    [serumdates, idx_sort] = sort(ratTable.SessionDate); %sort dates
    disp(serumdates)
    ratE2 = ratTable.EstradiolConc(idx_sort);
    E2_sess_overrats = [E2_sess_overrats; ratE2];
    stages_sess_overrats = [stages_sess_overrats; stages_rat];

end

%process
y = delta_sess_overrats;
x = E2_sess_overrats;
stages = stages_sess_overrats;

%exclude sessions where delta is NaN
y_good = y(~isnan(y));
x_good = x(~isnan(y));
stages_good = stages(~isnan(y));

%Binned (no estrus)
%exclude outliers
ub = mean(y_good, 'omitnan') + 1.5*std(y_good, 'omitnan');
lb = mean(y_good, 'omitnan') - 1.5*std(y_good, 'omitnan');
bad = find(y_good>=ub | y_good<=lb);
y_good(bad) = [];
x_good(bad) = [];
stages_good(bad) = [];

%bin data
xbin = linspace(min(x_good), max(x_good), 5);
ybin = NaN(1, length(xbin));
ybiner = NaN(1, length(xbin));
for j = 2:length(ybin)
    these = find(x_good>xbin(j-1) & x_good<=xbin(j));
    ybin(j) = mean(y_good(these));
    ybiner(j) = std(y_good(these))./sqrt(length(these));
end

end
