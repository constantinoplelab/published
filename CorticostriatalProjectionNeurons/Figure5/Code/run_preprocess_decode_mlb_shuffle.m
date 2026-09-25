
% creates the trainign and testing data for the projeciton neuron sthuffle 

function [] = run_preprocess_decode_mlb_shuffle(savename, usemlb, epoch, decodertype, usetorch)

    %location of analysis codebase
    rng(101)

    codepath = '~/projects/constantinoplelab/Analysis/';
    % add code paths
    addpath(genpath(codepath))
    addpath(genpath(codepath + "david"))
    addpath(genpath(codepath + "maggie"))
    
    if usetorch
        datapath = "/scratch/dh148/dynamics/data/maggie/";
    else
        datapath = "/Users/dhocker/projects/dynamics/data/maggie/";
    end
    b = load(datapath+"DS_projection_neurons_non-Stimulated.mat");
    
    %% decide on mlb or mlb2: update 20260107
    
    balance_blocks = decodertype;
    %% parse the data to get most likely block, time-averaged neurla responses on each trial
    output = dSTR_decode_parsedata(usetorch, epoch, balance_blocks);

    xvec = output.xvec;
    ntx = numel(xvec);
    ncond = 3;
    nfold = 20; %number of folds
    ntrials_percond = 20; %how many trials of each condition to withhold
    nsamps_perbucket_train = 500;
    nsamps_perbucket_test = 100; 
    
    %remove cells with only 8 high blok trials, which limit things
    cells2use_base = 1:60;
    cells2use_base(30:33) = [];
    
    cells2use = cells2use_base;
    ncell = numel(cells2use);
    
    
    %% build the testing dat
    disp('building testing data')
    trials_test_percell = nan(nfold, ncell, ncond, nsamps_perbucket_test);
    witheldtrials_test_percell = nan(nfold, ncell, ncond, ntrials_percond);
    % track the offers in each fold to see how balanced they are
    offers_test =  nan(nfold,ncell, ncond,nsamps_perbucket_test);
    
    % the testing data. do pseuodtrials for this too. 
    X_test = nan(nfold, ncell, ncond, nsamps_perbucket_test, ntx);
    Y_test = nan(nfold,ncond,nsamps_perbucket_test, ntx);
    
    counts_percat = nan(ncell,ncond);
    
    for m = 1:nfold
        for n = 1:ncell
            cellid = cells2use(n);
            conds_n = output.mostlikely_block{cellid};
            conds_n2 = output.secondlikely_block{cellid};
            if usemlb
                conds2use = conds_n;
            else
                conds2use = conds_n2;
            end
    
            rewardedtrials = b.all_S_DS{output.sess_ids(cellid)}.hits'==1;
    
            for k = 1:3
                trials_condk = find(conds2use == k & rewardedtrials); %ordered as high then low for 2nd most likely
                counts_percat(n,k) = numel(trials_condk);
                %disp([m,n,k,numel(trials_condk)])
                withheld_trials_kn = datasample(trials_condk, ntrials_percond, 'replace',false);          
                witheldtrials_test_percell(m,n,k,:) = withheld_trials_kn;
                trials_test_percell(m,n,k,:) = datasample(withheld_trials_kn,nsamps_perbucket_test,'replace',true);
                % track the offer size
                offers_test(m,n,k,:) = output.offer{cellid}(trials_test_percell(m,n,k,:));
    
                %populate testing data
                X_test(m,n,k,:,:) = output.X{cellid}(trials_test_percell(m,n,k,:),:);
                Y_test(m,k,:,:) = k;
            end
    
        end
    end
    
    
    %kill any nans. this will greatly slow the analysis since early and
    %late times will now be regressed
    X_test(isnan(X_test)) = 1e-8; % choose value that is identifiable
    
    %% build the training sets for each fold
    
    disp('building training data')
    ntx = numel(output.xvec); % number of regressions by timepoint
    X_training = nan(nfold,ncell, ncond,nsamps_perbucket_train,ntx); %training data
    Y_training = nan(nfold,ncond, nsamps_perbucket_train, ntx); % class
    trials_train_percell = nan(nfold,ncell,ncond,nsamps_perbucket_train);
    
    for m = 1:nfold
    
        for n = 1:ncell
            cellid = cells2use(n);
            %sess_idx_n = output.sess_ids(cellid);
    
            conds_n = output.mostlikely_block{cellid};
            conds_n2 = output.secondlikely_block{cellid};
            rewardedtrials = b.all_S_DS{output.sess_ids(cellid)}.hits'==1;
            if usemlb
                conds2use = conds_n;
            else
                conds2use = conds_n2;
            end
        
            % go through each category, remove samples from testing
            for k = 1:3
    
                mask = conds2use == k & rewardedtrials;
                allowed_cond_kn = find(mask);
                for kk = trials_test_percell(m,n,k,:)
                    allowed_cond_kn(allowed_cond_kn == kk) = [];
                end
    
                % sample with repalcement
                trials_train_percell(m,n,k,:) = datasample(allowed_cond_kn, nsamps_perbucket_train, 'replace',true);
                X_training(m,n,k,:,:) = output.X{cellid}(trials_train_percell(m,n,k,:),:);
                Y_training(m,k,:,:) = k;
    
            end
    
        end
    
    end
    %kill any nans. this will greatly slow the analysis since early and
    %late times will now be regressed
    X_training(isnan(X_training)) = 1e-8; % choose value that is identifiable
    
    %% save
    disp('saving')
    save(savename,'*');
    
    %savedir = '/Users/dhocker/projects/dynamics/data/maggie/';
    %if usemlb
    %    savename = strcat(savedir,'preprocess_mlb_quantile_shuffle_psthsubtract.mat');
    %else
    %    savename = strcat(savedir,'preprocess_mlb2_quantile_shuffle_psthsubtract.mat');
    %end
    

end