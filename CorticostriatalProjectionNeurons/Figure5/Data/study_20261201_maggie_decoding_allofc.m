% all OFC decoding. very long runtime

%location of analysis codebase
codepath = '~/projects/constantinoplelab/Analysis/';

% add code paths
addpath(genpath(codepath))
addpath(genpath(codepath + "david"))
addpath(genpath(codepath + "maggie"))
%% 
Etable = load('/Users/dhocker/projects/dynamics/data/maggie/EphysTable.mat')
ETable = Etable.ETable;

usemlb = true; %boolean to choose mlb or mlb2

%%
test = ETable{:,["recording_site"]};
mask1 = cellfun(@(x) strcmp(x,'OFC'),test,'UniformOutput',true);

test = ETable{:,["fiber_site"]};
mask2 = cellfun(@(x) strcmp(x,'DLS'),test,'UniformOutput',true);
%this does not matter

test = ETable{:,["protocol"]};
mask3 = cellfun(@(x) strcmp(x,'RWTautowait2'),test,'UniformOutput',true);
%mask3 = cellfun(@(x) strcmp(x,'RWTautowait2OptoTest'),test,'UniformOutput',true);
% does not currwently matter

test = ETable{:,["stimulation"]};
%mask4 = test==1; % stimulation
mask4 = test==0; % no stimulatio

savepath = ETable{:,["fullpath"]};
fname = ETable{:,["matfile"]};

% session number 1.

goodsessions = find(mask1 & mask2 & mask3 & mask4);
%goodsessions = find(mask1 & mask4);
nsess = numel(goodsessions);

disp('number of sessions')
disp(nsess)

%% 

loadmac = true;

fullnames = cell(nsess,1);

for j = 1:nsess

    if loadmac
        usedir = ETable{goodsessions(j),'savepath'}{1};
        usedir = replace(usedir,'\\constantinoplelab.cns.nyu.edu\server2\','/Volumes/server2/')
        usedir = replace(usedir,'\\constantinoplelab.cns.nyu.edu\server3\','/Volumes/server3/')
        usedir = replace(usedir,'\','/');

    else
        usedir = ETable{goodsessions(j),'savepath'}{1};
    end

    
    fullnames{j}  = strcat(usedir,ETable{goodsessions(j),'matfile'}{1});
    disp(ETable{goodsessions(j),'matfile'}{1})
end

% if .mat file was not present, remove from cell and get new nsess
goodsess = cellfun(@(x) numel(strfind(x,".mat")) > 0,fullnames);
fullnames = fullnames(goodsess);
nsess = numel(fullnames);

%%  get session data, count number neurons
Sall = cell(nsess,1);
epoch = 'reward';
%for j = 1:nsess
for j = 1:1
    disp(j)
    disp(fullnames{j})
    sj = struct();
    [output] = decode_parsedata_rawdat(fullnames{j},epoch);
    sj.output = output;
    sj.id = j;

    %check if there are enough samples per catetory: 
    ntrials_percond = 20; %number of trials per condition
    mlb_m = sum(sj.output.mostlikely_block == 1 & sj.output.rewarded_trials'==1);
    mlb_h = sum(sj.output.mostlikely_block == 2 & sj.output.rewarded_trials'==1);
    mlb_l = sum(sj.output.mostlikely_block == 3 & sj.output.rewarded_trials'==1);
    mlb2_m = sum(sj.output.secondlikely_block == 1 & sj.output.rewarded_trials'==1);
    mlb2_h = sum(sj.output.secondlikely_block == 2 & sj.output.rewarded_trials'==1);
    mlb2_l = sum(sj.output.secondlikely_block == 3 & sj.output.rewarded_trials'==1);

    if usemlb
        if mlb_m >= ntrials_percond && mlb_h >= ntrials_percond && mlb_l >= ntrials_percond
            Sall{j} = sj;
        else
            Sall{j} = [];
        end
    else
        if mlb2_m >= ntrials_percond && mlb2_h >= ntrials_percond && mlb2_l >= ntrials_percond
            Sall{j} = sj;
        else
            Sall{j} = [];
        end
    end
end

%remove sessions with less than 
badsess = cellfun(@(x) size(x,1),Sall,'UniformOutput',true) == 0;
Sall(badsess) = [];
nsess = numel(Sall);

%% choose sessions, parse raw data
epoch = 'reward';
%true params
ncell = 56; %number of projection neuron cells
nshuff = 20; %number of shuffles to do the decoding
nfold = 20;
rng(101);

%debug params
%ncell = 10; %number of projection neuron cells
%nshuff = 2; %number of shuffles to do the decoding
%nfold = 3;
rng(101);

% for every shuffle, use a different neuron set, which means
% also different session set
randsess_idx_flat = datasample(1:nsess,nshuff*ncell, 'replace',true);

    
%% get number of neurons per session, for every neuron across folds, sample

randsess_idx = reshape(randsess_idx_flat,nshuff, ncell);
%calculate num nuerons
nn = zeros(nsess,1);
for j = 1:nsess
    nn(j) = size(Sall{j}.output.X,1);
end

%get neuron for every shuffled
neuron_idx = zeros(nshuff, ncell);
for j = 1:nshuff
    for k = 1:ncell
        sess_id = randsess_idx(j,k);
        neuron_idx(j,k) = datasample(1:nn(sess_id),1);
    end
end

%% initialize save structures
xvec = Sall{1}.output.xvec;
ntx = numel(xvec);
ncond = 3;
ntrials_percond = 20; %how many trials of each condition to withhold
nsamps_perbucket_train = 500;
nsamps_perbucket_test = 100; 
dozscore = false;

preds_CV = nan(nshuff,nfold, ntx, ncond*nsamps_perbucket_test); %predictions on test data
preds_train = nan(nshuff,nfold, ntx, ncond*nsamps_perbucket_train); %predictions on training
true_CV = nan(nshuff,nfold, ntx,ncond*nsamps_perbucket_test); %true values, testing
true_trian = nan(nshuff,nfold, ntx,ncond*nsamps_perbucket_train); % true values, training
pvals = nan(nshuff,nfold, ntx,(ncond-1)*(ncell+1)); %coefficient pvalues, to help find good neurons
coeffs_model = nan(nshuff,nfold, ntx, (ncond-1)*(ncell+1));

coeftest_model = nan(nshuff,nfold,ntx);  % linear hypothesis test that all regressors are 0
Rsquared_model = nan(nshuff,nfold,ntx);  % r-squared of training data

%% for every shuffle, get a new  set of neurons, then run decoding
% LONG RUNTIME. saves at every new neuron set iteration, just in case of
% failure

for j_shuff = 1:nshuff

    disp("initializing data for this set of neurons: "+ num2str(j_shuff))
    % build cell arrays of MLB, MLB2, rewarded, hits for each sampled
    % neuron
    mlb_cell = cell(ncell,1);
    mlb2_cell = cell(ncell,1);
    offer_cell = cell(ncell,1);
    X_cell = cell(ncell,1);
    rewardedtrials_cell = cell(ncell,1);

    %get unique sessions for each neuron, populate data
    for j = 1:ncell 
        mlb_cell{j} = Sall{randsess_idx(j_shuff,j)}.output.mostlikely_block;
        mlb2_cell{j} = Sall{randsess_idx(j_shuff,j)}.output.secondlikely_block;
        offer_cell{j} = Sall{randsess_idx(j_shuff,j)}.output.offer;
        rewardedtrials_cell{j} = Sall{randsess_idx(j_shuff,j)}.output.rewarded_trials;
        
        %select the neuron for this shuffle, and this fold
        tmp = Sall{randsess_idx(j_shuff,j)}.output.X;  
        xvec = Sall{randsess_idx(j_shuff,j)}.output.xvec;
        xtmp = tmp{neuron_idx(j_shuff,j)};
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
            dat_mu = mean(mean(X_cell{cellid}(rewardedtrials,:),'omitnan'));
            dat_std = std(mean(X_cell{cellid}(rewardedtrials,:),2,'omitnan'));

            if usemlb
                conds2use = conds_n;
            else
                conds2use = conds_n2;
            end
    
            for k = 1:3
                disp([m,n,k])
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
                    dat_mu = 0;
                    dat_std = 1;
                    X_test(m,n,k,:,:) = X_cell{cellid}(trials_test_percell(m,n,k,:),:);
                end
    
                Y_test(m,k,:,:) = k;
            end
    
        end
    end

    % check testing data------------
    disp('check testing data. should be mlb2 = 3')
    test = squeeze(witheldtrials_test_percell(1,2,3,:));
    mlb_cell{2}(test)
    mlb2_cell{2}(test)

    %training data----------------------
    disp('building training data')
    X_training = nan(nfold,ncell, ncond,nsamps_perbucket_train,ntx); %training data
    Y_training = nan(nfold,ncond, nsamps_perbucket_train, ntx); % class
    trials_train_percell = nan(nfold,ncell,ncond,nsamps_perbucket_train);

    for m = 1:nfold
    
        for n = 1:ncell
            cellid = n;
            conds_n = mlb_cell{cellid};
            conds_n2 = mlb2_cell{cellid};
            rewardedtrials = rewardedtrials_cell{cellid}';
             % for z-scoring
            dat_mu = mean(mean(X_cell{cellid}(rewardedtrials,:),'omitnan'));
            dat_std = std(mean(X_cell{cellid}(rewardedtrials,:),2,'omitnan'));

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
                    dat_mu = 0;
                    dat_std = 1;
                    X_training(m,n,k,:,:) = X_cell{cellid}(trials_train_percell(m,n,k,:),:);
                end
                Y_training(m,k,:,:) = k;
    
            end
    
        end
    
    end

    % checkk training data
    disp('check training data')
    test = squeeze(trials_train_percell(1,2,3,:));
    mlb_cell{2}(test)
    mlb2_cell{2}(test)

    % check that data is structured as expected----------------------------
    disp('check regressor data')
    figure(188)
    clf
    hold on
    trial2use = trials_train_percell(1,2,2, 10); %fold 1, cell 2, high block, trial 10
    d = squeeze(X_training(1,2,2,10,:));
    sesstest = randsess_idx(j_shuff,2);
    celltest = neuron_idx(j_shuff,2);
    dtest = Sall{sesstest}.output.X{celltest}(trial2use,:);
    plot(d)
    plot(dtest+0.5)

    % do the regression----------------------------------------------------
    disp('doing regression')
    for m = 1:nfold
        for j = 1:ntx
            disp([j_shuff,m,j])
            Xtrain_mj = reshape(X_training(m,:,:,:,j),ncell,nsamps_perbucket_train*ncond)';
            Ytrain_mj = reshape(Y_training(m,:,:,j),1,nsamps_perbucket_train*ncond)';
    
            Xtest_mj = reshape(X_test(m,:,:,:,j),ncell,nsamps_perbucket_test*ncond)';
            Ytest_mj = reshape(Y_test(m,:,:,j),1,nsamps_perbucket_test*ncond)';
            try
                MnrMdl = fitmnr(Xtrain_mj,Ytrain_mj,'IterationLimit',100);
                coeftest_model(j_shuff,m,j) = MnrMdl.coefTest;
                Rsquared_model(j_shuff,m,j) = MnrMdl.Rsquared.Ordinary;
    
                preds_CV(j_shuff,m,j,:) = MnrMdl.predict(Xtest_mj);
                true_CV(j_shuff,m,j,:) = Ytest_mj;
    
                preds_train(j_shuff,m,j,:) = MnrMdl.predict(Xtrain_mj);
                true_trian(j_shuff,m,j,:) = Ytrain_mj;
    
                pvals(j_shuff,m,j,:) = MnrMdl.Coefficients{:,['pValue']};
                coeffs_model(j_shuff,m,j,:) = MnrMdl.Coefficients{:,['Value']};
                coeff_names = MnrMdl.CoefficientNames;
    
            catch
                %MnrMdl = fitmnr(Xtrain_mj,Ytrain_mj,'IterationLimit',100);
                disp('model bad. skipping')
            end
        end
    end


% end the shuffle here?
disp(datetime)
% save every shuff just in case
dosave = true;
savedir = '/Users/dhocker/projects/dynamics/results/maggie/';
if usemlb
    savename = strcat(savedir,'decode_mlb_alldat_allOFC_quantile.mat');
else
    savename = strcat(savedir,'decode_second_mlb_alldat_allOFC_quantile.mat');
end

if dosave
    save(savename,'X_training', 'X_test', 'Y_training', 'Y_test',...
        'preds_CV','true_CV', 'true_trian', 'coeff_names',...
        'pvals','coeftest_model','coeff_names','coeffs_model');
end

end


%% do some performance checks

%% plot average and sem predictive performance over time
idx_jshuff = 1;
preds_all = zeros(nshuff,nfold,ntx);

for n = 1:nshuff
    for m = 1:nfold
        for j = 1:ntx
            preds_all(n,m,j) = sum(preds_CV(n,m,j,:)==true_CV(n,m,j,:))/numel(true_CV(n,m,j,:));
        end
    end
end

accuracy_pred_means = squeeze(mean(preds_all,2,'omitnan'));
accuracy_pred_mean = mean(accuracy_pred_means,1,'omitnan');
accuracy_pred_sem = std(accuracy_pred_means,[],1,'omitnan')/sqrt(nshuff);

figure(353)
clf
hold on
shadedErrorBar(output.xvec,accuracy_pred_mean, accuracy_pred_sem)
hline(1/3,'k')
vline(0,'black')

title(strcat('prediction accuracy (2nd LB) on held out test data'))
xlabel('time from reward (s)')
ylabel('accuracy')
set(gca,'fontsize',15)


%% plot all traces individually
idx_jshuff = 1;
preds_all = zeros(nshuff,nfold,ntx);

for n = 1:nshuff
    for m = 1:nfold
        for j = 1:ntx
            preds_all(n,m,j) = sum(preds_CV(n,m,j,:)==true_CV(n,m,j,:))/numel(true_CV(n,m,j,:));
        end
    end
end

accuracy_pred_means = squeeze(mean(preds_all,2,'omitnan'));
accuracy_pred_mean = mean(accuracy_pred_means,1,'omitnan');
accuracy_pred_sem = std(accuracy_pred_means,[],1,'omitnan')/sqrt(nshuff);

figure(373)
clf
hold on

for j = 1:nshuff
plot(output.xvec, accuracy_pred_means(j,:),'linewidth',0.5,'color',[0.6,0.6,0.6])
end
plot(output.xvec,accuracy_pred_mean, 'linewidth',2,'color','k')
hline(1/3,'k')
vline(0,'black')

title(strcat('prediction accuracy (2nd LB) on held out test data'))
xlabel('time from reward (s)')
ylabel('accuracy')
set(gca,'fontsize',15)




