function EncodingWindowAnalysis = DARegressionAnalysis_CG(ratlist,...
    resample_freq, rebin_size, algn_vec, BestFits, Bstruct, Pstruct)

RPEs_over_rats = cell(length(ratlist), 1);
RPEs_last = cell(length(ratlist), 1);
lasttrialnums = -20:1:-1;
%get RPEs from photometry rats, last 20 trials
for rat = 1:length(ratlist)

    theseRPEs = BestFits.BestFit.(ratlist{rat}).All.RPE;
    thisratTrial = BestFits.BestFit.(ratlist{rat}).Train.ratTrial;
    RPEs_last{rat} = theseRPEs(ismember(thisratTrial.BlockPosition, lasttrialnums));
    RPEs_over_rats{rat} = theseRPEs;

end

for rat = 1:length(ratlist)

    %identify rat
    ratname = ratlist{rat};
    disp(['DA regression:' ratname])

    for aa  = 1:length(algn_vec)

        algn = algn_vec{aa};

        %subset to this event
        bdata = Bstruct.(ratname).(algn);
        pdata = Pstruct.(ratname).(algn);

        %get RPEs
        ratRPEs = RPEs_over_rats{rat};

        %Align dates of photometry and model
        ARat = BestFits.BestFit.(ratname).All.ratTrial;
        uniqueModelDates = ARat.date;
        mdldates_cell = cellstr(uniqueModelDates);
        mdldates = datetime(cellfun(@(sess) sess, mdldates_cell,...
            UniformOutput=false), InputFormat='dd-MMM-yyyy');
        rep_sess = find(diff(mdldates)==0);
        if ~isempty(rep_sess)
            uniqueModelDates(rep_sess, :) = [];
            ratRPEs(sum(ARat.ntrials(1:rep_sess-1))+1:sum(ARat.ntrials(1:rep_sess-1))+ARat.ntrials(rep_sess)) = [];
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
        RPEphot = ratRPEs(ismember(ModelDates_dt, sharedDates),:);
        MdlDatesphot = ModelDates_dt(ismember(ModelDates_dt, sharedDates),:);
        TrialNumsphot = ARat.trial_num(ismember(ModelDates_dt, sharedDates),:);
        BlockPosphot = ARat.BlockPosition(ismember(ModelDates_dt, sharedDates),:);
        Ntrialsphot = ARat.ntrials(ismember(datetime(uniqueModelDates), sharedDates),:);

        %subset
        %Model
        RPEs = RPEphot; %(logical(BLphot==1))
        MdlDates = MdlDatesphot; %(logical(BLphot==1))
        MdlTrialNums = TrialNumsphot; %(logical(BLphot==1))
        BlockPos = BlockPosphot; %(logical(BLphot==1))
        %Photometry
        bdata = bdataMdl; %((logical(bdataMdl.Block==1)),:)
        bdataDates = datetime(bdata.UniqueDay, 'InputFormat','uuuuMMdd');
        pdata = pdataMdl; %((logical(bdataMdl.Block==1)),:)

        %only keep mixed trials from model and photometry that match in trial number
        RPE_paireddown = NaN(size(bdata, 1), 1);
        BlockPos_paireddown = NaN(size(bdata, 1), 1);
        count = 0;
        for sess = 1:length(Ntrialsphot)
            %get trial numbers
            mdl_trialnums_sess = MdlTrialNums(ismember(MdlDates, sharedDates(sess)));
            phot_trialnums_sess = bdata.TrialNumber(ismember(...
                bdataDates, sharedDates(sess)));
            numtrials_phot = length(phot_trialnums_sess);

            %Get this session RPE
            RPE_thissess = RPEs(ismember(MdlDates, sharedDates(sess)));
            BlockPos_thissess = BlockPos(ismember(MdlDates, sharedDates(sess)));

            %keep RPE trials that are in photometry data
            mdltrials_in_phot = ismember(mdl_trialnums_sess, phot_trialnums_sess);
            RPE_paireddown(count+1:count+numtrials_phot) =...
                RPE_thissess(mdltrials_in_phot);
            BlockPos_paireddown(count+1:count+numtrials_phot) =...
                BlockPos_thissess(mdltrials_in_phot);

            phot_not_in_model = sum(~ismember(phot_trialnums_sess, mdl_trialnums_sess));
            if phot_not_in_model > 0
                sp(['photometry trials not in Mdl: '...
                    num2str(phot_not_in_model)])
            end

            count = count + length(phot_trialnums_sess);

        end

        %subset to last trials
        RPEs_thesetrials = RPE_paireddown(ismember(BlockPos_paireddown, lasttrialnums));
        da_mat = pdata(ismember(BlockPos_paireddown, lasttrialnums), :);
        bstruct_mat = bdata(ismember(BlockPos_paireddown, lasttrialnums), :);
        bstruct_mat.RPE = RPEs_thesetrials./max(abs(RPEs_thesetrials));

        [DAMat, times_resampled] = resample_da_data(da_mat, resample_freq);

        Regression = regress_da(DAMat, bstruct_mat); %times_resampled

        DA_Rebinned =...
            bin_da_by_var(DAMat, Regression, bstruct_mat, rebin_size);

        EncodingWindowAnalysis.(algn).DAMat = DAMat;
        EncodingWindowAnalysis.(algn).Times = times_resampled(1:end-1);
        EncodingWindowAnalysis.(algn).Regression = Regression;
        EncodingWindowAnalysis.(algn).RebinnedDA = DA_Rebinned;

    end


end

function [DAMat, times_resample] = resample_da_data(da_mat, resample_freq)
    times = linspace(-5, 10, size(da_mat,2));
    times_resample = -5:resample_freq:10;
    resample_indx = discretize(times, times_resample);

    DAMat = arrayfun(@(ii) mean(da_mat(:, resample_indx == ii), 2, 'omitnan'),...
        unique(resample_indx), 'UniformOutput', 0);
    DAMat = cell2mat(DAMat);
end

function EncodingWindow =...
        regress_da(da_mat_resamp, bstruct_all) %times_resampled

    EncodingWindow = struct;
    var = 'RPE';

    EncodingWindow.(var).betaMat = nan(2, size(da_mat_resamp, 2));
    EncodingWindow.(var).pvalMat = nan(2, size(da_mat_resamp, 2));

    for ii = 1:size(da_mat_resamp, 2)
        [EncodingWindow.(var).betaMat(:, ii), stats] =...
            robustfit(bstruct_all.(var), da_mat_resamp(:, ii));

        EncodingWindow.(var).pvalMat(:,ii) = stats.p;
    end

    % EncodingWindow.(var).betaMat(:, times_resampled>2.5) = nan;
    % EncodingWindow.(var).betaMat(:, times_resampled<0) = nan;

    [~, EncodingWindow.(var).I] =...
        max(abs(EncodingWindow.(var).betaMat(2,:)));
end

function [DA_Rebinned] = bin_da_by_var(da_mat, EncodingWindow,...
        bstruct_all, rebin_size)

    sem = @(mat) std(mat, [], 1, 'omitnan')./sqrt(size(mat, 1));

    DA_Rebinned = struct;
    DA_Rebinned.RebinSize = rebin_size;

    vars = fields(EncodingWindow);

    for v = 1:length(vars)

        var = vars{v};
        var_vec = bstruct_all.(var);

        binedges = -1:rebin_size:1;
        resample_indx = discretize(var_vec, binedges);

        mean_da = arrayfun(@(x) mean(da_mat(resample_indx==x, :), 1, 'omitnan'),...
            1:length(binedges)-1, 'UniformOutput', 0);

        sem_da = arrayfun(@(x) sem(da_mat(resample_indx==x, :)),...
            1:length(binedges)-1, 'UniformOutput', 0);

        DA_Rebinned.(var).DA_Mean = cell2mat(mean_da');
        DA_Rebinned.(var).DA_SEM =  cell2mat(sem_da');
        DA_Rebinned.(var).binedges = binedges;
    end


end

end