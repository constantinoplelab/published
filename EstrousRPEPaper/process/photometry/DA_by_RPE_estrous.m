function [pro_DA_binned, di_DA_binned, RPEbins_equallyspaced] =...
    DA_by_RPE_estrous(ratlist, Pstruct, Bstruct,...
    Proestrus_alphas, Diestrus_alphas, pro_ratList, di_ratList,...
    numbins, window, event)
%%% Dopamine AUC as a function of late alpha model-predicted RPE
% Just mixed blocks, last 20 trials, non-post-violations

%Get RPEs for photometry rats from mixed blocks and generate bins with
%equally distributed density
pro_RPEs = cell(length(ratlist), 1);
pro_RPEs_last = cell(length(ratlist), 1);
di_RPEs = cell(length(ratlist), 1);
di_RPEs_last = cell(length(ratlist), 1);
lasttrialnums = -20:1:-1;
%get RPEs from photometry rats
for rat = 1:length(ratlist)

    proratname = pro_ratList{contains(pro_ratList, ratlist{rat})};
    theseproRPEs = Proestrus_alphas.BestFit.(proratname).All.RPE;
    thisproratTrial = Proestrus_alphas.BestFit.(proratname).Train.ratTrial;
    pro_RPEs_last{rat} = theseproRPEs(ismember(thisproratTrial.BlockPosition, lasttrialnums));
    pro_RPEs{rat} = theseproRPEs;

    diratname = di_ratList{contains(di_ratList, ratlist{rat})};
    thesediRPEs = Diestrus_alphas.BestFit.(diratname).All.RPE;
    thisdiratTrial = Diestrus_alphas.BestFit.(diratname).Train.ratTrial;
    di_RPEs_last{rat} = thesediRPEs(ismember(thisdiratTrial.BlockPosition, lasttrialnums));
    di_RPEs{rat} = thesediRPEs;    

end
%Bin RPEs from last 20 trials
RPEs_rats = [pro_RPEs_last; di_RPEs_last];
[~, RPEbins_equallyspaced] = get_RPEbins(RPEs_rats, numbins);

%Get binned DA AUC at event, split by stage
pro_DA_binned = NaN(length(ratlist), numbins);
di_DA_binned = NaN(length(ratlist), numbins);
for rat = 1:length(ratlist)
    
    %identify rat
    ratname = ratlist{rat};
    disp(ratname)

    %load bdata and pdata
    bdata = Bstruct.(ratname).(event);
    pdata = Pstruct.(ratname).(event);

    %%PROESTRUS
    %get RPEs
    proratname = pro_ratList{contains(pro_ratList, ratname)};
    ratproRPEs = pro_RPEs{rat}; 

    %Align dates of photometry and model
    ARat = Proestrus_alphas.BestFit.(proratname).All.ratTrial;
    uniqueModelDates = ARat.date;
    mdldates_cell = cellstr(uniqueModelDates);
    mdldates = datetime(cellfun(@(sess) sess, mdldates_cell,...
        UniformOutput=false), InputFormat='dd-MMM-yyyy');
    %remove behavioral sessions from the same day (remove first instance of '09-Jul-2020' for G008 diestrus)
    rep_sess = find(diff(mdldates)==0);
    if ~isempty(rep_sess)
        uniqueModelDates(rep_sess, :) = [];
        pro_RPEs(sum(ARat.ntrials(1:rep_sess-1))+1:sum(ARat.ntrials(1:rep_sess-1))+ARat.ntrials(rep_sess)) = [];
        ARat.block(sum(ARat.ntrials(1:rep_sess-1))+1:sum(ARat.ntrials(1:rep_sess-1))+ARat.ntrials(rep_sess)) = [];
        ARat.trial_num(sum(ARat.ntrials(1:rep_sess-1))+1:sum(ARat.ntrials(1:rep_sess-1))+ARat.ntrials(rep_sess)) = [];
        ARat.ntrials(rep_sess) = [];
    end
    ModelDates = arrayfun(@(n)...
        repmat({uniqueModelDates(n, :)}, [1, ARat.ntrials(n)]),...
        1:length(ARat.ntrials), UniformOutput=false);
    ModelDates = [ModelDates{:}]';
    PDataDates = bdata.UniqueDay;
    ModelDates_dt = cellfun(@(x) datetime(x), ModelDates);
    PDataDates_dt = cellfun(@(x) datetime(x, 'InputFormat','uuuuMMdd'),...
        PDataDates);
    sharedDates = intersect(ModelDates_dt, PDataDates_dt);
    %photometry
    pdataMdl = pdata(ismember(PDataDates_dt, sharedDates),:);
    bdataMdl = bdata(ismember(PDataDates_dt, sharedDates),:);
    %model
    RPEphot = ratproRPEs(ismember(ModelDates_dt, sharedDates),:);
    MdlDatesphot = ModelDates_dt(ismember(ModelDates_dt, sharedDates),:);
    TrialNumsphot = ARat.trial_num(ismember(ModelDates_dt, sharedDates),:); 
    BlockPosphot = ARat.BlockPosition(ismember(ModelDates_dt, sharedDates),:);
    Ntrialsphot = ARat.ntrials(ismember(datetime(uniqueModelDates), sharedDates),:);

    %subset to mixed blocks
    %Model
    proRPEs_mixed = RPEphot; %(logical(BLphot==1))
    MdlDates_mixed = MdlDatesphot; %(logical(BLphot==1))
    MdlTrialNums_mixed = TrialNumsphot; %(logical(BLphot==1))
    BlockPos_mixed = BlockPosphot; %(logical(BLphot==1))
    %Photometry
    probdata_mixed = bdataMdl; %((logical(bdataMdl.Block==1)),:)
    bdataDates = datetime(probdata_mixed.UniqueDay, 'InputFormat','uuuuMMdd');
    bdataDates_mixed = bdataDates; %(logical(bdataMdl.Block==1))
    propdata_mixed = pdataMdl; %((logical(bdataMdl.Block==1)),:)

    %only keep mixed trials from model and photometry that match in trial number
    proRPE_paireddown = NaN(size(probdata_mixed, 1), 1);
    BlockPos_paireddown = NaN(size(probdata_mixed, 1), 1);
    count = 0;
    for sess = 1:length(Ntrialsphot)
        %get trial numbers
        mdl_trialnums_sess = MdlTrialNums_mixed(ismember(MdlDates_mixed, sharedDates(sess)));
        phot_trialnums_sess = probdata_mixed.TrialNumber(ismember(...
            bdataDates_mixed, sharedDates(sess)));
        numtrials_phot = length(phot_trialnums_sess);

        %Get this session RPE
        RPE_thissess = proRPEs_mixed(ismember(MdlDates_mixed, sharedDates(sess)));
        BlockPos_thissess = BlockPos_mixed(ismember(MdlDates_mixed, sharedDates(sess)));
        % Mdloptouts_thissess = Mdloptouts_mixed(ismember(MdlDates_mixed, sharedDates(sess)));

        %keep RPE trials that are in photometry data
        mdltrials_in_phot = ismember(mdl_trialnums_sess, phot_trialnums_sess);
        proRPE_paireddown(count+1:count+numtrials_phot) =...
            RPE_thissess(mdltrials_in_phot);
        BlockPos_paireddown(count+1:count+numtrials_phot) =...
            BlockPos_thissess(mdltrials_in_phot);        
     
        phot_not_in_model = sum(~ismember(phot_trialnums_sess, mdl_trialnums_sess));
        if phot_not_in_model > 0
            disp(['photometry trials not in Mdl: '...
                num2str(phot_not_in_model)])
        end

        count = count + length(phot_trialnums_sess);
        
    end

    %subset to last trials
    proRPEs_thesetrials = proRPE_paireddown(ismember(BlockPos_paireddown, lasttrialnums));
    propdata_thesetrials = propdata_mixed(ismember(BlockPos_paireddown, lasttrialnums), :);

    %%DIESTRUS
    %get RPEs
    diratname = di_ratList{contains(di_ratList, ratname)};
    ratdiRPEs = di_RPEs{rat};

    %Align dates of photometry and model
    ARat = Diestrus_alphas.BestFit.(diratname).All.ratTrial;
    uniqueModelDates = ARat.date;
    mdldates_cell = cellstr(uniqueModelDates);
    mdldates = datetime(cellfun(@(sess) sess, mdldates_cell,...
        UniformOutput=false), InputFormat='dd-MMM-yyyy');
    %remove behavioral sessions from the same day (remove first instance of '09-Jul-2020' for G008 diestrus)
    rep_sess = find(diff(mdldates)==0);
    if ~isempty(rep_sess)
        uniqueModelDates(rep_sess, :) = [];
        ratdiRPEs(sum(ARat.ntrials(1:rep_sess-1))+1:sum(ARat.ntrials(1:rep_sess-1))+ARat.ntrials(rep_sess)) = [];
        ARat.block(sum(ARat.ntrials(1:rep_sess-1))+1:sum(ARat.ntrials(1:rep_sess-1))+ARat.ntrials(rep_sess)) = [];
        ARat.trial_num(sum(ARat.ntrials(1:rep_sess-1))+1:sum(ARat.ntrials(1:rep_sess-1))+ARat.ntrials(rep_sess)) = [];
        ARat.ntrials(rep_sess) = [];
    end
    ModelDates = arrayfun(@(n)...
        repmat({uniqueModelDates(n, :)}, [1, ARat.ntrials(n)]),...
        1:length(ARat.ntrials), UniformOutput=false);
    ModelDates = [ModelDates{:}]';
    PDataDates = bdata.UniqueDay;
    ModelDates_dt = cellfun(@(x) datetime(x), ModelDates);
    PDataDates_dt = cellfun(@(x) datetime(x, 'InputFormat','uuuuMMdd'),...
        PDataDates);
    sharedDates = intersect(ModelDates_dt, PDataDates_dt);
    %photometry
    pdataMdl = pdata(ismember(PDataDates_dt, sharedDates),:);
    bdataMdl = bdata(ismember(PDataDates_dt, sharedDates),:);
    %model
    RPEphot = ratdiRPEs(ismember(ModelDates_dt, sharedDates),:);
    MdlDatesphot = ModelDates_dt(ismember(ModelDates_dt, sharedDates),:);
    TrialNumsphot = ARat.trial_num(ismember(ModelDates_dt, sharedDates),:); 
    BlockPosphot = ARat.BlockPosition(ismember(ModelDates_dt, sharedDates),:);
    Ntrialsphot = ARat.ntrials(ismember(datetime(uniqueModelDates), sharedDates),:);

    %subset to mixed blocks
    %Model
    diRPEs_mixed = RPEphot; %(logical(BLphot==1))
    MdlDates_mixed = MdlDatesphot; %(logical(BLphot==1))
    MdlTrialNums_mixed = TrialNumsphot; %(logical(BLphot==1))
    BlockPos_mixed = BlockPosphot; %(logical(BLphot==1))
    %Photometry
    dibdata_mixed = bdataMdl; %((logical(bdataMdl.Block==1)),:)
    bdataDates = datetime(dibdata_mixed.UniqueDay, 'InputFormat','uuuuMMdd');
    bdataDates_mixed = bdataDates; %(logical(dibdata_mixed.Block==1))
    dipdata_mixed = pdataMdl; %((logical(bdataMdl.Block==1)),:)

    %only keep mixed trials from model and photometry that match in trial number
    diRPE_paireddown = NaN(size(dibdata_mixed, 1), 1);
    BlockPos_paireddown = NaN(size(dibdata_mixed, 1), 1);
    count = 0;
    for sess = 1:length(Ntrialsphot)
        %get trial numbers
        mdl_trialnums_sess = MdlTrialNums_mixed(ismember(MdlDates_mixed, sharedDates(sess)));
        phot_trialnums_sess = dibdata_mixed.TrialNumber(ismember(...
            bdataDates_mixed, sharedDates(sess)));
        numtrials_phot = length(phot_trialnums_sess);

        %Get this session RPE
        RPE_thissess = diRPEs_mixed(ismember(MdlDates_mixed, sharedDates(sess)));
        BlockPos_thissess = BlockPos_mixed(ismember(MdlDates_mixed, sharedDates(sess)));        
        % Mdloptouts_thissess = Mdloptouts_mixed(ismember(MdlDates_mixed, sharedDates(sess)));

        %keep RPE trials that are in photometry data
        mdltrials_in_phot = ismember(mdl_trialnums_sess, phot_trialnums_sess);
        diRPE_paireddown(count+1:count+numtrials_phot) =...
            RPE_thissess(mdltrials_in_phot);
        BlockPos_paireddown(count+1:count+numtrials_phot) =...
            BlockPos_thissess(mdltrials_in_phot);   
        % Mdloptouts_paireddown(count+1:count+numtrials_phot) =...
        %     Mdloptouts_thissess(mdltrials_in_phot);

        phot_not_in_model = sum(~ismember(phot_trialnums_sess, mdl_trialnums_sess));
        if phot_not_in_model > 0
            disp(['photometry trials not in Mdl: '...
                num2str(phot_not_in_model)])
        end

        count = count + length(phot_trialnums_sess);
        
    end

    %subset to last trials
    diRPEs_thesetrials = diRPE_paireddown(ismember(BlockPos_paireddown, lasttrialnums));
    dipdata_thesetrials = dipdata_mixed(ismember(BlockPos_paireddown, lasttrialnums), :);

    %get times
    T = linspace(-5, 10, size(pdataMdl, 2)-1);
    [~ , tzero] = min(abs(T+0));
    [~ , tafter] = min(abs(T-window));

    %get mean DA by bin, split by stage
    for bin = 2:numbins+1

        theseproRPEs = proRPEs_thesetrials > RPEbins_equallyspaced(bin-1)... %RPEbins_density
            & proRPEs_thesetrials < RPEbins_equallyspaced(bin); %RPEbins_density
        if sum(theseproRPEs) > 1
            prodmat = mean(propdata_thesetrials(theseproRPEs, 2:end), 'omitnan');
            proAUCDA = trapz(T(:,tzero:tafter), prodmat(:,tzero:tafter));
            pro_DA_binned(rat, bin-1) = proAUCDA;
        else
            pro_DA_binned(rat, bin-1) = NaN;
        end

        thesediRPEs = diRPEs_thesetrials > RPEbins_equallyspaced(bin-1)... %RPEbins_density
            & diRPEs_thesetrials < RPEbins_equallyspaced(bin); %RPEbins_density
        if sum(thesediRPEs) > 1
            didmat = mean(dipdata_thesetrials(thesediRPEs, 2:end), 'omitnan');
            diAUCDA = trapz(T(:,tzero:tafter), didmat(:,tzero:tafter));
            di_DA_binned(rat, bin-1) = diAUCDA;
        else
            di_DA_binned(rat, bin-1) = NaN;
        end

    end

end


end

