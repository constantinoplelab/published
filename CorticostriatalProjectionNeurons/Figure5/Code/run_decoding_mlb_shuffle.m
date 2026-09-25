
function [] = run_decoding_mlb_shuffle(dataname, savename, shuff_idx)

    % for each shuffle, run code to decode
    rng(shuff_idx) % which shuffle to do

    %location of analysis codebase
    codepath = '~/projects/constantinoplelab/Analysis/';    
    % add code paths
    addpath(genpath(codepath))
    addpath(genpath(codepath + "david"))
    addpath(genpath(codepath + "maggie"))
    
    % set dimensionality things
    ntx = 61;
    ncond = 3;
    nfold = 20; %number of folds
    nsamps_perbucket_train = 500;
    nsamps_perbucket_test = 100; 


    f = load(dataname);
    X_training = f.X_training;
    Y_training = f.Y_training;
    X_test = f.X_test;
    Y_test = f.Y_test;

    ncell = size(X_test,2);
    
    % TODO: don't save the shuffle index
    preds_CV = nan(nfold, ntx, ncond*nsamps_perbucket_test); %predictions on test data
    preds_train = nan(nfold, ntx, ncond*nsamps_perbucket_train); %predictions on training
    true_CV = nan(nfold, ntx,ncond*nsamps_perbucket_test); %true values, testing
    true_trian = nan(nfold, ntx,ncond*nsamps_perbucket_train); % true values, training
    pvals = nan(nfold, ntx,(ncond-1)*(ncell+1)); %coefficient pvalues, to help find good neurons
    coeffs_model = nan(nfold, ntx, (ncond-1)*(ncell+1));
    
    coeftest_model = nan(nfold,ntx);  % linear hypothesis test that all regressors are 0
    Rsquared_model = nan(nfold,ntx);  % r-squared of training data


    % create shuffle X data
    Y_training_n = Y_training;
    % rather than sampling from class labels, just shuffle everything 
    shuff_labels_flat = datasample( [1,2,3], ncond*nsamps_perbucket_train);
    shuff_labels = reshape(shuff_labels_flat,ncond, nsamps_perbucket_train);
    
    % do the shuffle
    for j = 1:ncond
        for k = 1:nsamps_perbucket_train
            Y_training_n(:,j,k,:) = shuff_labels(j,k);
        end
    end

    % now iterate over the fold
    for m = 1:nfold
        for j = 1:ntx
            disp([m,j])
            Xtrain_mj = reshape(X_training(m,:,:,:,j),ncell,nsamps_perbucket_train*ncond)';
            Ytrain_mj = reshape(Y_training_n(m,:,:,j),1,nsamps_perbucket_train*ncond)';
    
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
        %end
    
    end

    %% save. name should include shuff_idx number
    
    save(savename,'preds_CV','true_CV', 'true_trian', 'coeff_names',...
            'pvals','coeftest_model','coeff_names','coeffs_model');

end