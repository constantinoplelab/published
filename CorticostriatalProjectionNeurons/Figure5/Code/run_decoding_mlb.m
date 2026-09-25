% most likely or 2nd likely block, using all data, not just mixed block as most likely
% this form will sutract the overall mean firing rate, conditioned on reward, from each neuron
% as a way of removing reward contributions, and potentially reducing
% accuracy of MLB2, which reviewers had concerns about


function [] = run_decoding_mlb(savename,usemlb, epoch, decodertype, usetorch)
    %% runs the mostlikely (and 2nd most likely) block decoder on projection neuron data
    % savename: (str) name of file to save results
    % uselmb: (bool) if true, do most likely block decoding. if false, 2nd
    %   most likely block
    % epoch: (str) task epoch. "reward","coff","son"
    % decodertype: (str) if "psth" will subtract off volume-specific psth
    %   from firing rates used to do regression. if "none", will use raw
    %   firing rates


    % set random seed
    rng(101)

    % the subtraction type. doing empirical mean before stratifying isn't
    % advised, so do the theoretical form of psth
    %balanceblocks = 'psth';

    % load code and data
    %location of analysis codebase
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
    
    %% parse the data to get most likely block, time-averaged neurla responses on each trial
    output = dSTR_decode_parsedata(usetorch, epoch, decodertype);
    
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
    trials_test_percell = nan(nfold, ncell, ncond, nsamps_perbucket_test);
    witheldtrials_test_percell = nan(nfold, ncell, ncond, ntrials_percond);
    
    % the testing data. do pseuodtrials for this too. 
    X_test = nan(nfold, ncell, ncond, nsamps_perbucket_test, ntx);
    Y_test = nan(nfold,ncond,nsamps_perbucket_test, ntx);
    
    counts_percat = nan(ncell,ncond);

    % track the offers in each fold to see how balanced they are
    offers_test =  nan(nfold,ncell, ncond,nsamps_perbucket_test);

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
                withheld_trials_kn = datasample(trials_condk, ntrials_percond, 'replace',false);          
                witheldtrials_test_percell(m,n,k,:) = withheld_trials_kn;
                trials_test_percell(m,n,k,:) = datasample(withheld_trials_kn,nsamps_perbucket_test,'replace',true);
    
                %populate testing data
                X_test(m,n,k,:,:) = output.X{cellid}(trials_test_percell(m,n,k,:),:);
                Y_test(m,k,:,:) = k+1;

                % track the offer size
                offers_test(m,n,k,:) = output.offer{cellid}(trials_test_percell(m,n,k,:));
            end
    
        end
    end
    %kill any nans. this will greatly slow the analysis since early and
    %late times will now be regressed
    X_test(isnan(X_test)) = 1e-8; % choose value that is identifiable

    %% build the training sets for each fold
  
    ntx = numel(output.xvec); % number of regressions by timepoint
    X_training = nan(nfold,ncell, ncond,nsamps_perbucket_train,ntx); %training data
    Y_training = nan(nfold,ncond, nsamps_perbucket_train, ntx); % class
    trials_train_percell = nan(nfold,ncell,ncond,nsamps_perbucket_train);
    offers_train =  nan(nfold,ncell, ncond,nsamps_perbucket_train);
    
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
                Y_training(m,k,:,:) = k+1;

                % track the offer size
                offers_train(m,n,k,:) = log2(output.offer{cellid}(trials_train_percell(m,n,k,:)));
    
            end
    
        end
    
    end
    
    %kill any nans. this will greatly slow the analysis since early and
    %late times will now be regressed
    X_training(isnan(X_training)) = 1e-8; % choose value that is identifiable
    
    
    %% for every fold, do nt regressions
    
    preds_CV = nan(nfold, ntx, ncond*nsamps_perbucket_test); %predictions on test data
    preds_train = nan(nfold, ntx, ncond*nsamps_perbucket_train); %predictions on training
    true_CV = nan(nfold, ntx,ncond*nsamps_perbucket_test); %true values, testing
    true_trian = nan(nfold, ntx,ncond*nsamps_perbucket_train); % true values, training
    pvals = nan(nfold, ntx,(ncond-1)*(ncell+1)); %coefficient pvalues, to help find good neurons
    coeffs_model = nan(nfold, ntx, (ncond-1)*(ncell+1));
    
    coeftest_model = nan(nfold,ntx);  % linear hypothesis test that all regressors are 0
    Rsquared_model = nan(nfold,ntx);  % r-squared of training data
    
    %for m = 1:1
    for m = 1:nfold
        for j = 1:ntx
            disp([m,j])
            Xtrain_mj = reshape(X_training(m,:,:,:,j),ncell,nsamps_perbucket_train*ncond)';
            Ytrain_mj = reshape(Y_training(m,:,:,j),1,nsamps_perbucket_train*ncond)';
    
            Xtest_mj = reshape(X_test(m,:,:,:,j),ncell,nsamps_perbucket_test*ncond)';
            Ytest_mj = reshape(Y_test(m,:,:,j),1,nsamps_perbucket_test*ncond)';


            try
                MnrMdl = fitmnr(Xtrain_mj,Ytrain_mj,'IterationLimit',100);
                coeftest_model(m,j) = MnrMdl.coefTest;
                Rsquared_model(m,j) = MnrMdl.Rsquared.Ordinary;
    
                preds_CV(m,j,:) = MnrMdl.predict(Xtest_mj);
                true_CV(m,j,:) = Ytest_mj;
    
                preds_train(m,j,:) = MnrMdl.predict(Xtrain_mj);
                true_trian(m,j,:) = Ytrain_mj;
    
                pvals(m,j,:) = MnrMdl.Coefficients{:,'pValue'};
                coeffs_model(m,j,:) = MnrMdl.Coefficients{:,'Value'};
                coeff_names = MnrMdl.CoefficientNames;
    
            catch
                %MnrMdl = fitmnr(Xtrain_mj,Ytrain_mj,'IterationLimit',100);
                disp('model bad. skipping')
            end
    
            
    
        end
    end
    
    %% save

    save(savename,'X_training', 'X_test', 'Y_training', 'Y_test',...
        'preds_CV','true_CV', 'true_trian', 'coeff_names',...
        'pvals','coeftest_model','coeff_names','coeffs_model', ...
        'offers_test','offers_train',...
        'trials_train_percell','trials_test_percell');

end





