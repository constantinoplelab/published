% figure-generating code for Maggie's paper. includes the following
% most likely block and 2nd most likely block decoding

%% load code path and data.

%location of analysis codebase
codepath = '~/projects/constantinoplelab/Analysis/';

% add code paths
addpath(genpath(codepath))
addpath(genpath(codepath + "david"))
addpath(genpath(codepath + "maggie"))

%for saving pds correctly so illustrator can modify them
set(0, 'DefaultFigureRenderer', 'painters');

% DS projection neuron data location
datadir = "Z:\david\maggie\";
savedir = "Z:\david\maggie\";

fname = datadir + "DS_projection_neurons_non-Stimulated.mat";
b = load(fname);

% parse the DS projections to get beliefs
output = dSTR_decode_parsedata();

%% load relevant pre-processed data: the regression results for each case

dprimetype = 'all'
usemlb = false; % if true, most likely block. if false, 2nd most likely

if usemlb
    savename_shuff = strcat(savedir,'decode_mlb_alldat_dprime_shuffle_',dprimetype,'_quantile.mat');
    savename = strcat(savedir,'decode_mlb_alldat_dprime_',dprimetype,'_quantile.mat');
    savename_ofc = strcat(savedir,'decode_mlb_alldat_allOFC_quantile.mat');
    savename_table = 'Z:\david\maggie\MLB.xlsx';
else
    savename_shuff = strcat(savedir,'decode_second_mlb_alldat_dprime_shuffle_',dprimetype,'_quantile.mat');
    savename = strcat(savedir,'decode_second_mlb_alldat_dprime_',dprimetype,'_quantile.mat');
    savename_ofc = strcat(savedir,'decode_second_mlb_alldat_allOFC_quantile.mat');
    savename_table = 'Z:\david\maggie\2MLB.xlsx';
end

d_shuff = load(savename_shuff);
d = load(savename);
d_ofc = load(savename_ofc);

%% plot average and sem predictive performance over time
nshuff = 20;
nfold = 20;
ntx = 61;

xvec_ofc = -4:0.2:4;
ntx_ofc = numel(xvec_ofc);

%calculate the shuffle mean and sem
preds_all = zeros(nshuff,nfold,ntx);
for n = 1:nshuff
    for m = 1:nfold
        for j = 1:ntx
            preds_all(n,m,j) = sum(d_shuff.preds_CV(n,m,j,:)==d_shuff.true_CV(n,m,j,:))/numel(d_shuff.true_CV(n,m,j,:));
        end
    end
end

accuracy_pred_means = squeeze(mean(preds_all,2,'omitnan'));
accuracy_shuff_mean = mean(accuracy_pred_means,1,'omitnan');
accuracy_shuff_sem = 1.96*std(accuracy_pred_means,[],1,'omitnan');


%calculate the ofc result
preds_all_ofc = zeros(nshuff,nfold,ntx_ofc);
for n = 1:nshuff
    for m = 1:nfold
        for j = 1:ntx_ofc
            preds_all_ofc(n,m,j) = sum(d_ofc.preds_CV(n,m,j,:)==d_ofc.true_CV(n,m,j,:))/numel(d_ofc.true_CV(n,m,j,:));
        end
    end
end

accuracy_pred_means_ofc = squeeze(mean(preds_all_ofc,2,'omitnan'));
accuracy_ofc_mean = mean(accuracy_pred_means_ofc,1,'omitnan');
%accuracy_ofc_mean = median(accuracy_pred_means_ofc,1,'omitnan');
accuracy_ofc_sem = std(accuracy_pred_means_ofc,[],1,'omitnan')/sqrt(20); %number of shuffles

% the real data
preds_all = zeros(nfold,ntx);
for m = 1:nfold
    for j = 1:ntx
        preds_all(m,j) = sum(d.preds_CV(m,j,:)==d.true_CV(m,j,:))/numel(d.true_CV(m,j,:));
    end
end

accuracy_pred_mean = mean(preds_all,1,'omitnan');
accuracy_pred_sem = std(preds_all,[],1,'omitnan')/sqrt(nfold);



figure(353)
clf
hold on
shadedErrorBar(output.xvec,accuracy_pred_mean, accuracy_pred_sem,'lineprops',{'color','red'})
shadedErrorBar(xvec_ofc,accuracy_ofc_mean, accuracy_ofc_sem, 'lineprops',{'color','blue'})
shadedErrorBar(output.xvec,accuracy_shuff_mean, accuracy_shuff_sem, 'lineprops',{'color',[0.05,0.05,0.05]})

if usemlb
    title(strcat(dprimetype,', prediction accuracy (MLB) on held out test data'))
else
    title(strcat(dprimetype,', prediction accuracy (2nd LB) on held out test data'))
end

xlabel('time from reward (s)')
ylabel('accuracy')

set(gca,'fontsize',15)
xlim([-1,2.5])

% do the significance tests comparing all OFC vs. dSTR
xtest = -1:0.2:2.5;
ntest = numel(xtest);
pvals = zeros(1,ntest);
for j = 1:numel(xtest)

    [~,idx_ofc] = min(abs(xvec_ofc-xtest(j)));
    [~,idx] = min(abs(output.xvec-xtest(j)));
    dat_ofc = accuracy_pred_means_ofc(:,idx_ofc);
    dat = preds_all(:,idx);
    
    pvals(j) = signrank(dat_ofc,dat);
end
pval2 = pvals;
pvals(pvals > 0.05) = nan;
pgood = find(pvals < 0.05);
plot(xtest(pgood),0.8*ones(size(pgood)),'*', 'linewidth', 0.5, 'color','k')
ylim([0,1.0])

legend('true','all ofc','shuffle', 'significant')
hline(1/3,'k--')
vline(0,'k--')

times = -1:0.2:2.5;
T = table(times', pval2','VariableNames',{'Time','p value'});
writetable(T,savename_table)
