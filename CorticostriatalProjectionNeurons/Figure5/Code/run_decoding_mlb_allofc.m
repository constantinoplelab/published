% all OFC decoding. very long runtime. does conditional psth subtraction


function [] = run_decoding_mlb_allofc(dataname,savename, usemlb, shuff2use)

    %location of analysis codebase
    codepath = '~/projects/constantinoplelab/Analysis/';
    rng(101)
    
    % add code paths
    addpath(genpath(codepath))
    addpath(genpath(codepath + "david"))
    addpath(genpath(codepath + "maggie"))

    % load data explicitly
    f = load(dataname);
    Sall = f.Sall;
    disp(savename)
    
    %% choose sessions, parse raw data
    %true params
    nsess = numel(Sall);
    ncell = 56; %number of projection neuron cells
    %nshuff = 20; %number of shuffles to do the decoding
    nshuff = 100;
    nfold = 20;
    
    %debug params
    %ncell = 10; %number of projection neuron cells
    %nshuff = 2; %number of shuffles to do the decoding
    %nfold = 3;
    rng(101);
    
    % some dimesionality things
    xvec = Sall{1}.output.xvec;
    ntx = numel(xvec);
    ncond = 3;
    ntrials_percond = 20; %how many trials of each condition to withhold
    nsamps_perbucket_train = 500;
    nsamps_perbucket_test = 100; 
    dozscore = false;


    %% get number of neurons per session, for every neuron across folds, sample
    %make sure the sessions have enough trials per condition. remove any that do not
    ntrials_percond_mlb = zeros(nsess,3);
    ntrials_percond_mlb2 = zeros(nsess,3);

    for j = 1:nsess
        for k = 1:3
            ntrials_percond_mlb(j,k) = sum(Sall{j}.output.rewarded_trials == 1 & Sall{j}.output.mostlikely_block' == k);
            ntrials_percond_mlb2(j,k) = sum(Sall{j}.output.rewarded_trials == 1 & Sall{j}.output.secondlikely_block' == k);
        end
    end
    
    
    % sessions where there are 20 or less trials per conditions
    sess2keep_mlb = sum(ntrials_percond_mlb > ntrials_percond,2) == 3 & sum(ntrials_percond_mlb2 > ntrials_percond,2) == 3;
    Sall = Sall(sess2keep_mlb);

    nsess = numel(Sall);
    disp('number of sessions')
    disp(nsess)
    
    
    %for every shuffle, use a different neuron set, which means
    % also different session set
    randsess_idx_flat = datasample(1:nsess,nshuff*ncell, 'replace',true);
        
    randsess_idx = reshape(randsess_idx_flat,nshuff, ncell);
    %calculate num nuerons
    nn = zeros(nsess,1);
   
    for j = 1:nsess
        nn(j) = size(Sall{j}.output.X,1);
    end
    
    %get neuron for every shuffled. avoid 1hz neurons
    neuron_idx = zeros(nshuff, ncell);
    for j = 1:nshuff
        for k = 1:ncell
            sess_id = randsess_idx(j,k);
            foundcell = false;
            while foundcell == false
                candidate_neuron = datasample(1:nn(sess_id),1);
                if Sall{sess_id}.output.meanrates(candidate_neuron) > 1.0
                    neuron_idx(j,k) = candidate_neuron;
                    foundcell = true;
                else
                    disp('cell with fr < 1hz suggested. trying again')
                    disp(Sall{sess_id}.output.meanrates(candidate_neuron))
                end
            end

            %neuron_idx(j,k) = datasample(1:nn(sess_id),1);
        end
    end
    
    %% initialize save structures. fix this to only include a single entry,
    % not the full, inflated nshuff one
    
    % TODO: don't save entire shuffle
    preds_CV = nan(nfold, ntx, ncond*nsamps_perbucket_test); %predictions on test data
    preds_train = nan(nfold, ntx, ncond*nsamps_perbucket_train); %predictions on training
    true_CV = nan(nfold, ntx,ncond*nsamps_perbucket_test); %true values, testing
    true_trian = nan(nfold, ntx,ncond*nsamps_perbucket_train); % true values, training
    pvals = nan(nfold, ntx,(ncond-1)*(ncell+1)); %coefficient pvalues, to help find good neurons
    coeffs_model = nan(nfold, ntx, (ncond-1)*(ncell+1));
    
    coeftest_model = nan(nfold,ntx);  % linear hypothesis test that all regressors are 0
    Rsquared_model = nan(nfold,ntx);  % r-squared of training data
    % TODO: what else to add
    offers_test =  nan(nfold,ncell, ncond,nsamps_perbucket_test);
    
    %% redo the full k-fold with random draws of neurons from ofc. single run
    % LONG RUNTIME. saves at every new neuron set iteration, just in case of
    % failure
    %for samp_idx = 1:nshuff
    for samp_idx = shuff2use

        disp(strcat('beginning ofc sample: ',num2str(samp_idx)));
        txt2add = strcat('_',num2str(samp_idx),'.mat');
        savename_shuff = replace(savename,'.mat',txt2add);
    
        % build cell arrays of MLB, MLB2, rewarded, hits for each sampled
        % neuron
        mlb_cell = cell(ncell,1);
        mlb2_cell = cell(ncell,1);
        offer_cell = cell(ncell,1);
        X_cell = cell(ncell,1);
        rewardedtrials_cell = cell(ncell,1);
        offers_cell = cell(ncell,1);
    
        %get unique sessions for each neuron, populate data
        for j = 1:ncell 
            mlb_cell{j} = Sall{randsess_idx(samp_idx,j)}.output.mostlikely_block;
            mlb2_cell{j} = Sall{randsess_idx(samp_idx,j)}.output.secondlikely_block;
            offer_cell{j} = Sall{randsess_idx(samp_idx,j)}.output.offer;
            rewardedtrials_cell{j} = Sall{randsess_idx(samp_idx,j)}.output.rewarded_trials;
            offers_cell{j} = Sall{randsess_idx(samp_idx,j)}.output.rewarded_trials;
            
            %select the neuron for this shuffle, and this fold
            tmp = Sall{randsess_idx(samp_idx,j)}.output.X;  
            %xvec = Sall{randsess_idx(samp_idx,j)}.output.xvec;
            xtmp = tmp{neuron_idx(samp_idx,j)};
            X_cell{j} = xtmp;
        end
    
        %testing data------------------------------
        disp('building testing data')
        trials_test_percell = nan(nfold, ncell, ncond, nsamps_perbucket_test);
        witheldtrials_test_percell = nan(nfold, ncell, ncond, ntrials_percond);
        
        % the testing data. do pseuodtrials for this too. 
        X_test = nan(nfold, ncell, ncond, nsamps_perbucket_test, ntx);
        Y_test = nan(nfold,ncond,nsamps_perbucket_test, ntx);
        
        counts_percat = nan(ncell,ncond);
        
        for m = 1:nfold
            for n = 1:ncell
        
                % get likely blocks and rewarded trials
                cellid = n;
                conds_n = mlb_cell{cellid};
                conds_n2 = mlb2_cell{cellid};
                rewardedtrials = rewardedtrials_cell{cellid}';
                % for z-scoring
                if dozscore
                    dat_mu = mean(mean(X_cell{cellid}(rewardedtrials,:),'omitnan'));
                    dat_std = std(mean(X_cell{cellid}(rewardedtrials,:),2,'omitnan'));
                end
    
                if usemlb
                    conds2use = conds_n;
                else
                    conds2use = conds_n2;
                end
        
                for k = 1:3
                    trials_condk = find(conds2use == k & rewardedtrials); %ordered as high then low for 2nd most likely
                    counts_percat(n,k) = numel(trials_condk);
                    %disp([m,n,k,numel(trials_condk)])
                    withheld_trials_kn = datasample(trials_condk, ntrials_percond, 'replace',false);          
                    witheldtrials_test_percell(m,n,k,:) = withheld_trials_kn;
                    trials_test_percell(m,n,k,:) = datasample(withheld_trials_kn,nsamps_perbucket_test,'replace',true);
        
                    %populate testing data. z-score
                    if dozscore
                        X_test(m,n,k,:,:) = (X_cell{cellid}(trials_test_percell(m,n,k,:),:)-dat_mu)/dat_std;
                    else
                        X_test(m,n,k,:,:) = X_cell{cellid}(trials_test_percell(m,n,k,:),:);
                    end
        
                    Y_test(m,k,:,:) = k;

                    % track the offer size
                    offers_test(m,n,k,:) = offers_cell{cellid}(trials_test_percell(m,n,k,:));
                end
        
            end
        end
    
    
        %training data----------------------
        disp('building training data')
        X_training = nan(nfold,ncell, ncond,nsamps_perbucket_train,ntx); %training data
        Y_training = nan(nfold,ncond, nsamps_perbucket_train, ntx); % class
        trials_train_percell = nan(nfold,ncell,ncond,nsamps_perbucket_train);
        offers_train =  nan(nfold,ncell, ncond,nsamps_perbucket_train);
    
        for m = 1:nfold
        
            for n = 1:ncell
                cellid = n;
                conds_n = mlb_cell{cellid};
                conds_n2 = mlb2_cell{cellid};
                rewardedtrials = rewardedtrials_cell{cellid}';
                 % for z-scoring
                if dozscore
                    dat_mu = mean(mean(X_cell{cellid}(rewardedtrials,:),'omitnan'));
                    dat_std = std(mean(X_cell{cellid}(rewardedtrials,:),2,'omitnan'));
                end
    
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
                    if dozscore
                        X_training(m,n,k,:,:) = (X_cell{cellid}(trials_train_percell(m,n,k,:),:)-dat_mu)/dat_std;
                    else
                        X_training(m,n,k,:,:) = X_cell{cellid}(trials_train_percell(m,n,k,:),:);
                    end
                    Y_training(m,k,:,:) = k;
        
                end

                % track the offer size
                offers_train(m,n,k,:) = offers_cell{cellid}(trials_train_percell(m,n,k,:));
        
            end   
        end
    

        %kill any nans. this will greatly slow the analysis since early and
        %late times will now be regressed
        X_training(isnan(X_training)) = 1e-8; % choose value that is identifiable

    
        % do the regression----------------------------------------------------
        disp('doing regression')
        for m = 1:nfold
            for j = 1:ntx
                disp([samp_idx,m,j])
               
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
    
    
    % end the shuffle here?
    disp(datetime)
    save(savename_shuff,...
        'preds_CV','true_CV', 'true_trian', 'coeff_names',...
        'pvals','coeftest_model','coeff_names','coeffs_model', ...
        'offers_test','offers_train',...
        'trials_train_percell','trials_test_percell');

    
    end
end


