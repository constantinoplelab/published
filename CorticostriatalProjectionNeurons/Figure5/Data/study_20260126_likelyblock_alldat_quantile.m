% most likely or 2nd likely block, using all data, not just mixed block as most likely

%% load code and data
%location of analysis codebase
codepath = '~/projects/constantinoplelab/Analysis/';

% add code paths
addpath(genpath(codepath))
addpath(genpath(codepath + "david"))
addpath(genpath(codepath + "maggie"))

datapath = "/Users/dhocker/projects/dynamics/data/maggie/";
b = load(datapath+"DS_projection_neurons_non-Stimulated.mat");
c  = load(datapath + "DS_projection_neurons_non-Stimulated2.mat");

%% decide if most likely or 2nd most likely 
usemlb = false;

%% parse the data to get most likely block, time-averaged neurla responses on each trial
output = dSTR_decode_parsedata();

X_percell = output.X;
mostlikelyblock = output.mostlikely_block;
secondmostlikelyblock = output.secondlikely_block;
alltrials = output.alltrials;

sess_ids = output.sess_ids;
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

tidx = 26; % t = 0.4 s after reward
dprimes2compare = c.signed_dprime_DS.Rew(cells2use_base,tidx);

%% subselect cells by d prime, but use top/bottom percentile to match numbers

figure(38)
clf
hold on
histogram(dprimes2compare,linspace(-1,1,21))
xlabel('d prime')
ylabel('counts')
title(strcat('d prime values of opto neurons at t = ',num2str(xvec(tidx)),'s'))
set(gca,'fontsize',15)

%clip some of the d prime values close to zero?
boundary = 0.2; % top and bottom 20%

[~,sortidx] = sort(dprimes2compare);
nct = numel(dprimes2compare); %total number of cells
ncell = floor(nct*0.2); % number of cells to use


dprimetype = 'high';

switch dprimetype
    case 'low'
        cells2use = cells2use_base(sortidx(1:ncell));
        ncell = numel(cells2use);
    case 'high'
        cells2use = cells2use_base(sortidx(nct-ncell+1:nct));
        ncell = numel(cells2use);

    case 'null'
        cells2use = cells2use_base(sortidx(ncell+1:nct-ncell));
        ncell = numel(cells2use);

    case 'all'
        cells2use = cells2use_base;
        ncell = numel(cells2use);
end


%% build the testing dat
trials_test_percell = nan(nfold, ncell, ncond, nsamps_perbucket_test);
witheldtrials_test_percell = nan(nfold, ncell, ncond, ntrials_percond);

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

            %populate testing data
            X_test(m,n,k,:,:) = output.X{cellid}(trials_test_percell(m,n,k,:),:);
            Y_test(m,k,:,:) = k+1;
        end

    end
end

%% check that trials are correct condition

% should be 1, then 3
test = squeeze(witheldtrials_test_percell(1,2,3,:));
output.mostlikely_block{cells2use(2)}(test)
output.secondlikely_block{cells2use(2)}(test)

%% build the training sets for each fold


ntx = numel(output.xvec); % number of regressions by timepoint
X_training = nan(nfold,ncell, ncond,nsamps_perbucket_train,ntx); %training data
Y_training = nan(nfold,ncond, nsamps_perbucket_train, ntx); % class
trials_train_percell = nan(nfold,ncell,ncond,nsamps_perbucket_train);

for m = 1:nfold

    for n = 1:ncell
        cellid = cells2use(n);
        sess_idx_n = output.sess_ids(cellid);

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

        end

    end

end

%% check trials
% should be 1, then 2
test = squeeze(trials_train_percell(1,2,2,:));
output.mostlikely_block{cells2use(2)}(test)
output.secondlikely_block{cells2use(2)}(test)

%% check the pca of a given training and testing fold, colored by mlb

fold = 1;
tidx = 22; %

dat_train_high = X_training(fold,:,1,:,tidx);
dat_train_low = X_training(fold,:,2,:,tidx);
dat_train_allcond = reshape(X_training(fold,:,:,:,tidx),ncell,ncond*nsamps_perbucket_train);

disp('here')

dat_test_high = X_test(fold,:,1,:,tidx);
dat_test_low = X_test(fold,:,2,:,tidx);
dat_test_allcond = reshape(X_test(fold,:,:,:,tidx),ncell,ncond*nsamps_perbucket_test);

[coeff,score,latent,tsquared,explained,mu] = pca(dat_train_allcond');

figure(291)
clf
plot(cumsum(explained),'linewidth',2)
xlabel('PC')
ylabel('cumulative variance explained')
title('cumulative variance explained training data')
set(gca,'fontsize',15)
hline(90,'k')



[coeff_test,score_test,latent_test,tsquared_test,explained_test,mu_test] = pca(dat_test_allcond');

figure(303)
clf
plot(cumsum(explained_test),'linewidth',2)
xlabel('PC')
ylabel('cumulative variance explained')
title('cumulative variance explained testing data')
set(gca,'fontsize',15)
hline(90,'k')



%%

figure(310)
clf
hold on

PC_train_high = squeeze(dat_train_high)'*coeff(:,1:2);
PC_train_low = squeeze(dat_train_low)'*coeff(:,1:2);

scatter(PC_train_high(:,1),PC_train_high(:,2), 'markerfacecolor','red', 'markeredgecolor','red')
scatter(PC_train_low(:,1),PC_train_low(:,2), 'markerfacecolor','blue', 'markeredgecolor','blue')
xlabel('pc 1')
ylabel('pc 2')
title('PC projections of training data, colored by 2nd most likely block')
set(gca,'fontsize',15)


figure(330)
clf
hold on

%PC_test_mixed = squeeze(dat_test_mixed)'*coeff_test(:,1:2);
%PC_test_high = squeeze(dat_test_high)'*coeff_test(:,1:2);
%PC_test_low = squeeze(dat_test_low)'*coeff_test(:,1:2);

PC_test_high = squeeze(dat_test_high)'*coeff(:,1:2);
PC_test_low = squeeze(dat_test_low)'*coeff(:,1:2);

scatter(PC_test_high(:,1),PC_test_high(:,2), 'markerfacecolor','red', 'markeredgecolor','red')
scatter(PC_test_low(:,1),PC_test_low(:,2), 'markerfacecolor','blue', 'markeredgecolor','blue')
xlabel('pc 1')
ylabel('pc 2')
title('PC projections of test data, colored by 2nd most likely block')
set(gca,'fontsize',15)


%% check that data is structured as expected, 
% and that reshaping preserves structure

figure(188)
clf
hold on
trial2use = trials_train_percell(1,2,2, 10); %fold 1, cell 20, high block, trial 10
d = squeeze(X_training(1,2,2,10,:));
dtest = b.all_SU_DLS_200{cells2use(2)}.hmat.Rew(trial2use,:);
plot(d)
plot(dtest+0.5)


%% for every fold, do nt regressions

preds_CV = nan(nfold, ntx, ncond*nsamps_perbucket_test); %predictions on test data
preds_train = nan(nfold, ntx, ncond*nsamps_perbucket_train); %predictions on training
true_CV = nan(nfold, ntx,ncond*nsamps_perbucket_test); %true values, testing
true_trian = nan(nfold, ntx,ncond*nsamps_perbucket_train); % true values, training
pvals = nan(nfold, ntx,(ncond-1)*(ncell+1)); %coefficient pvalues, to help find good neurons
coeffs_model = nan(nfold, ntx, (ncond-1)*(ncell+1));

coeftest_model = nan(nfold,ntx);  % linear hypothesis test that all regressors are 0
Rsquared_model = nan(nfold,ntx);  % r-squared of training data


%TODO: will need to get params
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
dosave = true;
doload = false;

savedir = '/Users/dhocker/projects/dynamics/results/maggie/';
if usemlb
    savename = strcat(savedir,'decode_mlb_alldat_dprime_',dprimetype,'_quantile.mat');
else
    savename = strcat(savedir,'decode_second_mlb_alldat_dprime_',dprimetype,'_quantile.mat');
end

if dosave
    save(savename,'X_training', 'X_test', 'Y_training', 'Y_test',...
        'preds_CV','true_CV', 'true_trian', 'coeff_names',...
        'pvals','coeftest_model','coeff_names','coeffs_model');
end

if doload
    load(savename)
end

%% look at performance fold 1

preds_test = squeeze(preds_CV(1,:,:));
true_test = squeeze(true_CV(1,:,:));

preds_accuracy = sum(preds_test == true_test,2)/(3*nsamps_perbucket_test);



figure(353)
clf
hold on
plot(output.xvec,preds_accuracy)

title('prediction accuracy on held out test data pseudotrials, single fold')
xlabel('time from reward (s)')
ylabel('accuracy')
set(gca,'fontsize',15)

%% confusion matrix

fold = 4;
t_idx = 26;

disp(sum(preds_CV(fold,t_idx,:)==true_CV(fold,t_idx,:))/numel(true_CV(fold,t_idx,:)))

C = confusionmat(squeeze(true_CV(fold,t_idx,:)), squeeze(preds_CV(fold,t_idx,:)));
figure(328)
clf
hold off


txtcell = {'mixed','high','low'};

confusionchart(C,txtcell);
set(gca,'fontsize',15)
title(strcat('confusion matrix mlb: fold=',num2str(fold),', time from reward = ',num2str(output.xvec(t_idx))))




%% plot average and sem predictive performance over time

preds_all = zeros(nfold,ntx);
for m = 1:nfold
    for j = 1:ntx
        preds_all(m,j) = sum(preds_CV(m,j,:)==true_CV(m,j,:))/numel(true_CV(m,j,:));
    end
end

accuracy_pred_mean = mean(preds_all,1,'omitnan');
accuracy_pred_sem = std(preds_all,[],1,'omitnan')/sqrt(nfold);

figure(353)
clf
hold on
shadedErrorBar(output.xvec,accuracy_pred_mean, accuracy_pred_sem)
hline(1/3,'k')
vline(0,'black')

title(strcat(dprimetype,', prediction accuracy (2nd LB) on held out test data'))
xlabel('time from reward (s)')
ylabel('accuracy')
set(gca,'fontsize',15)

%% same, for training data

preds_all = zeros(nfold,ntx);
for m = 1:nfold
    for j = 1:ntx
        preds_all(m,j) = sum(preds_train(m,j,:)==true_trian(m,j,:))/numel(true_trian(m,j,:));
    end
end

accuracy_pred_mean = mean(preds_all,1,'omitnan');
accuracy_pred_sem = std(preds_all,[],1,'omitnan');

figure(387)
clf
hold on
shadedErrorBar(output.xvec,accuracy_pred_mean, accuracy_pred_sem)

title('prediction accuracy on training pseudotrials')
xlabel('time from reward (s)')
ylabel('accuracy')
set(gca,'fontsize',15)
hline(1/3)
vline(0.5)


%% look at coefficients and their p values to find neurons that care
t_idx = 28;
accuracy_pval_mean = median(log10(pvals(:,t_idx,:)),1,'omitnan');
accuracy_pval_sem = std(log10(pvals(:,t_idx,:)),[],1,'omitnan');

figure(348)
clf
hold on
shadedErrorBar(1:ncell+1,accuracy_pval_mean, accuracy_pval_sem )

title(strcat('median +sem log pvals at t = ',num2str(output.xvec(t_idx))))
xlabel('regressor')
ylabel('log(p)')
set(gca,'fontsize',15)
hline(log10(0.05),'k')


%% look at psth by volume for given neuron

cell_idx = 9;
dat = output.X{cells2use(cell_idx)};
dat_s = b.all_S_DS{b.all_index_DS{cells2use(cell_idx),4}};

colors = {[0.5,0,0.8],'blue',[0,0.7,0],[0.8,0.8,0],'red'};

vmask_5 = dat_s.RewardAmount == 5 & dat_s.hits==1;
vmask_10 = dat_s.RewardAmount == 10 & dat_s.hits==1;
vmask_20 = dat_s.RewardAmount == 20 & dat_s.hits==1;
vmask_40 = dat_s.RewardAmount == 40 & dat_s.hits==1;
vmask_80 = dat_s.RewardAmount == 80 & dat_s.hits==1;

figure(415)
clf
hold on
plot(xvec,mean(dat(vmask_5,:),1,'omitnan'),'linewidth',2,'color',colors{1})
plot(xvec,mean(dat(vmask_10,:),1,'omitnan'),'linewidth',2,'color',colors{2})
plot(xvec,mean(dat(vmask_20,:),1,'omitnan'),'linewidth',2,'color',colors{3})
plot(xvec,mean(dat(vmask_40,:),1,'omitnan'),'linewidth',2,'color',colors{4})
plot(xvec,mean(dat(vmask_80,:),1,'omitnan'),'linewidth',2,'color',colors{5})
vline(xvec(22),'r-.')
vline(0,'k-')

legend('5','10','20','40','80')
title(strcat('neuron ',num2str(cells2use(cell_idx)),', dprime type ',dprimetype))
xlabel('time to reward (s)')
ylabel('psth')
set(gca,'fontsize',15)

%% split the neurons by class and sign


names2use = coeff_names(good_reg);
nregressor = size(names2use,1);
ismixedvar = nan(nregressor,1);
reg_sign = nan(nregressor,1);
regnum = nan(nregressor,1);
for j = 1:nregressor
    ismixedvar(j) = numel(strfind(names2use(j,:),'_1')) > 0;
    reg_sign(j) = median(coeffs_model(:,t_idx,good_reg(j)),'omitnan') > 0;
    name_j = names2use{j};
    us_idx = strfind(name_j,'_');
    regnum(j) = str2num(name_j(2:us_idx-1));
end

% mixed neurons with positive weights
mixed_pos = ismixedvar & reg_sign;
mixed_neg = ismixedvar & ~reg_sign;

high_pos = ~ismixedvar & reg_sign;
high_neg =  ~ismixedvar & ~reg_sign;

% grab the cells that have mixed-voting and positive sign

neurons_mixedpos = cells2use(regnum(find(mixed_pos)));
didsp(neurons_mixedpos)

figure(422)
clf
hold on

neurons_mixed = cells2use(regnum(find(mixed_pos | mixed_neg)));
nn = numel(neurons_mixed);
for jj = 1:nn
    subplot(1,nn,jj)
    j = neurons_mixed(jj);
    % find trials to average
    mlb_j = output.mostlikely_block{j} == 1 & b.all_S_DS{output.sess_ids(j)}.hits' == 1;
    mlb_notj = output.mostlikely_block{j} ~= 1 & b.all_S_DS{output.sess_ids(j)}.hits' == 1;

    hmat_j = output.X{j};
    
    mean_j = mean(hmat_j(mlb_j,:),1,'omitnan');   
    sem_j = std(hmat_j(mlb_j,:),[],1,'omitnan')/sqrt(sum(mlb_j));

    mean_j_non = mean(hmat_j(mlb_notj,:),1,'omitnan');
    sem_j_non = std(hmat_j(mlb_notj,:),[],1,'omitnan')/sqrt(sum(mlb_notj));

    shadedErrorBar(output.xvec, mean_j,sem_j,'lineprops',{'color','k'});
    shadedErrorBar(output.xvec, mean_j_non,sem_j_non,'lineprops',{'color','b'});
    vline(0,'k')
    vline(output.xvec(t_idx),'k--')
    xlabel('time to reward (s)')
    ylabel('firing rate (z-scored)')
    title(strcat('neuron ',num2str(j),', MLB. mixed'))
    set(gca,'fontsize',15)
end


%% neuron 25 should be better for a high block than a low block. check it
figure(451)
clf
hold on

j = 25; %neuron
% find trials to average
mlb_j = output.mostlikely_block{j} == 1 & b.all_S_DS{output.sess_ids(j)}.hits' == 1;
mlb_notj = output.mostlikely_block{j} ~= 1 & b.all_S_DS{output.sess_ids(j)}.hits' == 1;

hmat_j = output.X{j};

mean_j = mean(hmat_j(mlb_j,:),1,'omitnan');   
sem_j = std(hmat_j(mlb_j,:),[],1,'omitnan')/sqrt(sum(mlb_j));

mean_j_non = mean(hmat_j(mlb_notj,:),1,'omitnan');
sem_j_non = std(hmat_j(mlb_notj,:),[],1,'omitnan')/sqrt(sum(mlb_notj));

shadedErrorBar(output.xvec, mean_j,sem_j,'lineprops',{'color','k'});
shadedErrorBar(output.xvec, mean_j_non,sem_j_non,'lineprops',{'color','b'});
vline(0,'k')
vline(output.xvec(t_idx),'k--')
xlabel('time to reward (s)')
ylabel('firing rate (z-scored)')
title(strcat('neuron ',num2str(j),', MLB. mixed'))
set(gca,'fontsize',15)

 %% neuron 2 should be different in mixed than low
 figure(451)
clf
hold on

    j = 2; %neuron
    % find trials to average
    mlb_j = output.mostlikely_block{j} == 2 & b.all_S_DS{output.sess_ids(j)}.hits' == 1;
    mlb_notj = output.mostlikely_block{j} ~= 2 & b.all_S_DS{output.sess_ids(j)}.hits' == 1;

    hmat_j = output.X{j};
    
    mean_j = mean(hmat_j(mlb_j,:),1,'omitnan');   
    sem_j = std(hmat_j(mlb_j,:),[],1,'omitnan')/sqrt(sum(mlb_j));

    mean_j_non = mean(hmat_j(mlb_notj,:),1,'omitnan');
    sem_j_non = std(hmat_j(mlb_notj,:),[],1,'omitnan')/sqrt(sum(mlb_notj));

    shadedErrorBar(output.xvec, mean_j,sem_j,'lineprops',{'color','r'});
    shadedErrorBar(output.xvec, mean_j_non,sem_j_non,'lineprops',{'color','b'});
    vline(0,'k')
    vline(output.xvec(t_idx),'k--')
    xlabel('time to reward (s)')
    ylabel('firing rate (z-scored)')
    title(strcat('neuron ',num2str(j),', MLB. mixed'))
    set(gca,'fontsize',15)



%% look at rsquared values for training data

rsquard_pred_mean = mean(Rsquared_model,1,'omitnan');
rsquared_pred_sem = std(Rsquared_model,[],1,'omitnan');

figure(402)
clf
hold on
shadedErrorBar(output.xvec,rsquard_pred_mean, rsquared_pred_sem)

title('r squared on traiing pseudotrials')
xlabel('time from reward (s)')
ylabel('R squared')
set(gca,'fontsize',15)

%% look at the coefficient test



coefftest_pred_mean = mean(log10(coeftest_model),1,'omitnan');
coefftest_pred_sem = std(log10(coeftest_model),[],1,'omitnan');

figure(419)
clf
hold on
shadedErrorBar(output.xvec,coefftest_pred_mean, coefftest_pred_sem)
hline(log10(0.05),'k')
title('log of coefficient test squared on traiing pseudotrials')
xlabel('time from reward (s)')
ylabel('R squared')
set(gca,'fontsize',15)


%% do LOOCV

ntrial = size(X_flat,1);
pred_LOOCV = nan(ntrial,1);
ncoeff = 122; %number of coefficients
pvals = nan(ntrial, ncoeff); % pvalues
coeffvals = nan(ntrial, ncoeff); 

for m = 1:ntrial
    disp(m)
    X_test = X_flat(m,:);
    Y_test = Y_flat(m,:);
    train_mask = true(ntrial,1);
    train_mask(m) = false;
    X_train = X_flat(train_mask, :);
    Y_train = Y_flat(train_mask);

    % multinomial logistic regression
    MnrMdl = fitmnr(X_train,Y_train);
    Ypred = predict(MnrMdl,X_test);
    
    pred_LOOCV(m) = Ypred;
    coeff = MnrMdl.Coefficients;
    coeffnames = MnrMdl.CoefficientNames;
    pvals(m,:) = coeff{:,'pValue'};
    coeffvals(m,:) = coeff{:,'Value'};
end

%% percent correct and confusion matrix?
disp(sum(pred_LOOCV==Y_flat)/numel(Y_flat))

C = confusionmat(Y_flat,pred_LOOCV);
figure(221)
clf
hold off

txtcell = {'mixed','high','low'};

confusionchart(C,txtcell)
set(gca,'fontsize',15)
title('confusion matrix for predicting most likely block. LOOCV')

%% find what is significant

% look at log pvalues
figure(161)
clf
hold on
plot(log10(pvals))
xlabel('coefficients')
ylabel('log p value')
hline(log10(0.05),'k')
title('p values of coefficients form all fits')
set(gca,'fontsize',15)




