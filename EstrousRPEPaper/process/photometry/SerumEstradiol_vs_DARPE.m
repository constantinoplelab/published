function [delta_good, E2_good, stages] =...
    SerumEstradiol_vs_DARPE(SerumTable, PstructSerum, BstructSerum)
% correlate estradiol expression to DA AUC at largest reward offer cue per session (putative largest RPE)

%set variables for for loop
ratlist = fields(BstructSerum);
deltaAUC_sess_overrats = [];
E2_sess_overrats = [];
stages_sess_overrats = [];
window = 0.5; %time after CPIn to be included in the area under the curve (AUC) calculation (0.5 s)
T = linspace(-5, 10, size(PstructSerum.('G127').('CPIn'), 2)-1); %time vector, just chose one of the rats' datasets here to create this from
[~ , tzero] = min(abs(T+0)); %which time bin contains the CPIn event
[~ , tafter] = min(abs(T-window)); %which time bin contains the end of the AUC window
[~ , tstartbl] = min(abs(T+0.05)); %for baseline correction below, starts 0.05 s before CPIn
for rat = 1:length(ratlist)

    ratname = ratlist{rat};
    disp(ratname)

    %get rat's photometry data
    pstruct = PstructSerum.(ratname).('CPIn');
    bdata = BstructSerum.(ratname).('CPIn');

    %max normalize based on the mean response in mixed blocks.
    mblocks = bdata.Block==1 & bdata.PrevTrialType~=2;
    mixed_resp = mean(pstruct(mblocks,:), 'omitnan');
    pstruct = pstruct./max(mixed_resp); %divide by max of CPIn response in mixed blocks

    %subset to serum sample sessions
    serum_sessiondates = SerumTable.SessionDate(strcmp(SerumTable.Rat, ratname));
    serum_sessiondates_unique = unique(serum_sessiondates); %orders them as well
    DAhsess_withserum_idx = cellfun(@(x)...
        ismember(datetime(x, 'InputFormat','uuuuMMdd'),...
        datetime(serum_sessiondates_unique)), bdata.UniqueDay);
    DAserumdates = unique(bdata.UniqueDay(DAhsess_withserum_idx));
    disp([num2str(length(serum_sessiondates_unique)) ' sessions captured'])
    disp(DAserumdates)
    %exclude noisy G138
    if strcmp(ratname, 'G138')
        lowSNR = ismember(DAserumdates, {'20240627', '20240701',...
            '20240709'});
        DAserumdates(lowSNR) = [];
    end

    %get DA AUC for RPE-like events
    for sess = 1:length(DAserumdates)

        high_trials = bdata.Block==2 & bdata.PrevTrialType~=2 & ...
            strcmp(bdata.UniqueDay, DAserumdates(sess));
        pdata_high = pstruct(high_trials, 2:end);
        high_mat = mean(pdata_high, 'omitnan');
        baseline = min(high_mat(:,tstartbl:tzero));
        high_bc = high_mat - baseline;
        AUC_high = trapz(T(:,tzero:tafter), high_bc(:,tzero:tafter));

        low_trials = bdata.Block==3 & bdata.PrevTrialType~=2 &...
            strcmp(bdata.UniqueDay, DAserumdates(sess));
        pdata_low = pstruct(low_trials, 2:end);
        low_mat = mean(pdata_low, 'omitnan');
        baseline = min(low_mat(:,tstartbl:tzero));
        low_bc = low_mat - baseline;
        AUC_low = trapz(T(:,tzero:tafter), low_bc(:,tzero:tafter));

        deltaAUC_sess = AUC_high - AUC_low; 
    
        deltaAUC_sess_overrats = [deltaAUC_sess_overrats; deltaAUC_sess];
        ratE2 = SerumTable.EstradiolConc(strcmp(SerumTable.Rat, ratname) &...
            ismember(datetime(SerumTable.SessionDate),...
            datetime(DAserumdates{sess}, 'InputFormat','uuuuMMdd')), :); 
        E2_sess_overrats = [E2_sess_overrats; ratE2];
        stages_sess_overrats = [stages_sess_overrats;...
            unique(bdata.stage(strcmp(bdata.UniqueDay,...
            DAserumdates(sess))))];

    end

end

delta = deltaAUC_sess_overrats;
E2 = E2_sess_overrats;
stages = stages_sess_overrats;

% %remove estrus
% delta(strcmp(stages, 'Estrus')) = nan;

%y(isoutlier(y)) = nan;
delta_good = delta(~isnan(delta)); %exclude sessions where AUC is NaN
E2_good = E2(~isnan(delta)); %exclude sessions where AUC is NaN
stages = stages(~isnan(delta));

end