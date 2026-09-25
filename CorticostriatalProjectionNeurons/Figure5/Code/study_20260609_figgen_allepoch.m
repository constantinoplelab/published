% check the dSTR for all epochs, and also compare decoder types

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
% TODO: MAGGIE change these paths to point to server or your local machine
% if you want to run this code
datadir = '\\constantinoplelab.cns.nyu.edu\server\david\maggie\data\'; % "/Users/dhocker/projects/dynamics/data/maggie/";
savedir = 'Z:\david\maggie\data\'; %"/Users/dhocker/projects/dynamics/results/maggie/";

fname = datadir + "DS_projection_neurons_non-Stimulated.mat";
b = load(fname);

% parse the DS projections to get beliefs
usetorch = false;
output = dSTR_decode_parsedata(usetorch);
dprimetype = 'all';


% xvec = output.xvec;
output.xvec = -4:0.2:8;
xvec_ofc = -4:0.2:4;
ntx = numel(xvec_ofc);
ntx_ofc = numel(xvec_ofc);
nfold = 20;

%% set the decodertype and stratify type

decodertype = 'psth2';
dct = decodertype; %easier to type 
is_stratified = false;


%% load relevant pre-processed data: the regression results for each case

if is_stratified

    s_reward = strcat(savedir,'20260604_stratified/mlb_reward_',dct,'.mat');
    s_reward_ofc = strcat(savedir,'allofc_stratified/postprocess_mlb_allofc_reward_',dct,'_stratified.mat');
    s_reward_shuffle = strcat(savedir,'shuffle_stratified/postprocess_mlb_shuffle_reward_',dct,'_stratified.mat');

    s_coff = strcat(savedir,'20260604_stratified/mlb_coff_',dct,'.mat');
    s_coff_ofc = strcat(savedir,'allofc_stratified/postprocess_mlb_allofc_coff_',dct,'_stratified.mat');
    s_coff_shuffle = strcat(savedir,'shuffle_stratified/postprocess_mlb_shuffle_coff_',dct,'_stratified.mat');

    s_son = strcat(savedir,'20260604_stratified/mlb_son_none.mat');
    s_son_ofc = strcat(savedir,'allofc_stratified/postprocess_mlb_allofc_son_',dct,'_stratified.mat');
    s_son_shuffle = strcat(savedir,'shuffle_stratified/postprocess_mlb_shuffle_son_',dct,'_stratified.mat');


    % mlb2
    s_reward2 = strcat(savedir,'20260604_stratified/mlb2_reward_',dct,'.mat');
    s_reward_ofc2 = strcat(savedir,'allofc_stratified/postprocess_mlb2_allofc_reward_',dct,'_stratified.mat');
    s_reward_shuffle2 = strcat(savedir,'shuffle_stratified/postprocess_mlb2_shuffle_reward_',dct,'_stratified.mat');

    s_coff2 = strcat(savedir,'20260604_stratified/mlb2_coff_',dct,'.mat');
    s_coff_ofc2 = strcat(savedir,'allofc_stratified/postprocess_mlb2_allofc_coff_',dct,'_stratified.mat');
    s_coff_shuffle2 = strcat(savedir,'shuffle_stratified/postprocess_mlb2_shuffle_coff_',dct,'_stratified.mat');

    s_son2 = strcat(savedir,'20260604_stratified/mlb2_son_',dct,'.mat');
    s_son_ofc2 = strcat(savedir,'allofc_stratified/postprocess_mlb2_allofc_son_',dct,'_stratified.mat');
    s_son_shuffle2 = strcat(savedir,'shuffle_stratified/postprocess_mlb2_shuffle_son_',dct,'_stratified.mat');

else
    s_reward = strcat(savedir,'20260604_decode_projection/mlb_reward_',dct,'.mat');
    s_reward_ofc = strcat(savedir,'allofc/postprocess_mlb_allofc_reward_',dct,'.mat');
    s_reward_shuffle = strcat(savedir,'shuffle/postprocess_mlb_shuffle_reward_',dct,'.mat');

    s_coff = strcat(savedir,'20260604_decode_projection/mlb_coff_',dct,'.mat');
    s_coff_ofc = strcat(savedir,'allofc/postprocess_mlb_allofc_coff_',dct,'.mat');
    s_coff_shuffle = strcat(savedir,'shuffle/postprocess_mlb_shuffle_coff_',dct,'.mat');

    s_son = strcat(savedir,'20260604_decode_projection/mlb_son_',dct,'.mat');
    s_son_ofc = strcat(savedir,'allofc/postprocess_mlb_allofc_son_',dct,'.mat');
    s_son_shuffle = strcat(savedir,'shuffle/postprocess_mlb_shuffle_son_',dct,'.mat');

    % mlb 2
    s_reward2 = strcat(savedir,'20260604_decode_projection/mlb2_reward_',dct,'.mat');
    s_reward_ofc2 = strcat(savedir,'allofc/postprocess_mlb2_allofc_reward_',dct,'.mat');
    s_reward_shuffle2 = strcat(savedir,'shuffle/postprocess_mlb2_shuffle_reward_',dct,'.mat');

    s_coff2 = strcat(savedir,'20260604_decode_projection/mlb2_coff_',dct,'.mat');
    s_coff_ofc2 = strcat(savedir,'allofc/postprocess_mlb2_allofc_coff_',dct,'.mat');
    s_coff_shuffle2 = strcat(savedir,'shuffle/postprocess_mlb2_shuffle_coff_',dct,'.mat');

    s_son2 = strcat(savedir,'20260604_decode_projection/mlb2_son_',dct,'.mat');
    s_son_ofc2 = strcat(savedir,'allofc/postprocess_mlb2_allofc_son_',dct,'.mat');
    s_son_shuffle2 = strcat(savedir,'shuffle/postprocess_mlb2_shuffle_son_',dct,'.mat');
   

end
    

% load the data

d_reward = load(s_reward);
d_reward_ofc = load(s_reward_ofc);
d_reward_shuffle = load(s_reward_shuffle); 
d_coff = load(s_coff);
d_coff_ofc = load(s_coff_ofc);
d_coff_shuffle = load(s_coff_shuffle);
d_son = load(s_son);
d_son_ofc = load(s_son_ofc); % data not present yet
d_son_shuffle = load(s_son_shuffle); %data not present yet

d_reward2 = load(s_reward2);
d_reward_ofc2 = load(s_reward_ofc2);
d_reward_shuffle2 = load(s_reward_shuffle2); 
d_coff2 = load(s_coff2);
d_coff_ofc2 = load(s_coff_ofc2);
d_coff_shuffle2 = load(s_coff_shuffle2);
d_son2 = load(s_son2);
d_son_ofc2 = load(s_son_ofc2); % data not present yet
d_son_shuffle2 = load(s_son_shuffle2); %data not present yet

%% set some testing parameters
xl_reward = [-1,3];
xl_son = [-1,3];
xl_coff = [0,3];
alpha = 0.05;

%% calcualte the accuracy for the projection neurons

cfun = @(x) causal_smooth(x,5);
%mfun = @(x) mean(x,1,'omitnan');
mfun = @(x) median(x,1,'omitnan');
%cfun = @(x) x; %no smoothing

epoch = 'reward';
d = d_reward;
preds_all = zeros(nfold,ntx);
for m = 1:nfold
    for j = 1:ntx
        preds_all(m,j) = sum(d.preds_CV(m,j,:)==d.true_CV(m,j,:))/numel(d.true_CV(m,j,:));
    end
end
accuracy_reward_psth_mean = cfun(mfun(preds_all));
accuracy_reward_psth_sem = cfun(std(preds_all,[],1,'omitnan')/sqrt(nfold));
%tests
d_sh = d_reward_shuffle;
d_ofc = d_reward_ofc;
accuracy_pred_means_ofc = mean(d_ofc.preds_all_ofc,2,'omitnan');
accuracy_pred_means_shuffle = mean(d_sh.preds_all_shuffle,2,'omitnan');
[pgood_mlb_reward, pv_reward_shuff, pv_reward_ofc] = sigtests(preds_all,accuracy_pred_means_ofc,accuracy_pred_means_shuffle,output.xvec, alpha);



d = d_reward2;
preds_all = zeros(nfold,ntx);
for m = 1:nfold
    for j = 1:ntx
        preds_all(m,j) = sum(d.preds_CV(m,j,:)==d.true_CV(m,j,:))/numel(d.true_CV(m,j,:));
    end
end
accuracy_reward_psth_mean2 = cfun(mfun(preds_all));
accuracy_reward_psth_sem2 = cfun(std(preds_all,[],1,'omitnan')/sqrt(nfold));
%tests
d_sh = d_reward_shuffle2;
d_ofc = d_reward_ofc2;
accuracy_pred_means_ofc = mean(d_ofc.preds_all_ofc,2,'omitnan');
accuracy_pred_means_shuffle = mean(d_sh.preds_all_shuffle,2,'omitnan');
[pgood_mlb2_reward, pv_reward_shuff2, pv_reward_ofc2] = sigtests(preds_all,accuracy_pred_means_ofc,accuracy_pred_means_shuffle,output.xvec, alpha);




% coff
d = d_coff;
preds_all = zeros(nfold,ntx);
for m = 1:nfold
    for j = 1:ntx
        preds_all(m,j) = sum(d.preds_CV(m,j,:)==d.true_CV(m,j,:))/numel(d.true_CV(m,j,:));
    end
end
accuracy_coff_psth_mean = cfun(mfun(preds_all));
accuracy_coff_psth_sem = cfun(std(preds_all,[],1,'omitnan')/sqrt(nfold));
%tests
d_sh = d_coff_shuffle;
d_ofc = d_coff_ofc;
accuracy_pred_means_ofc = mean(d_ofc.preds_all_ofc,2,'omitnan');
accuracy_pred_means_shuffle = mean(d_sh.preds_all_shuffle,2,'omitnan');
[pgood_mlb_coff, pv_coff_shuff, pv_coff_ofc] = sigtests(preds_all,accuracy_pred_means_ofc,accuracy_pred_means_shuffle,output.xvec, alpha);




d = d_coff2;
preds_all = zeros(nfold,ntx);
for m = 1:nfold
    for j = 1:ntx
        preds_all(m,j) = sum(d.preds_CV(m,j,:)==d.true_CV(m,j,:))/numel(d.true_CV(m,j,:));
    end
end
accuracy_coff_psth_mean2 = cfun(mfun(preds_all));
accuracy_coff_psth_sem2 = cfun(std(preds_all,[],1,'omitnan')/sqrt(nfold));
%tests
d_sh = d_coff_shuffle;
d_ofc = d_coff_ofc;
accuracy_pred_means_ofc = mean(d_ofc.preds_all_ofc,2,'omitnan');
accuracy_pred_means_shuffle = mean(d_sh.preds_all_shuffle,2,'omitnan');
[pgood_mlb2_coff, pv_coff_shuff2, pv_coff_ofc2] = sigtests(preds_all,accuracy_pred_means_ofc,accuracy_pred_means_shuffle,output.xvec, alpha);



% son
d = d_son;
preds_all = zeros(nfold,ntx);
for m = 1:nfold
    for j = 1:ntx
        preds_all(m,j) = sum(d.preds_CV(m,j,:)==d.true_CV(m,j,:))/numel(d.true_CV(m,j,:));
    end
end
accuracy_son_psth_mean = cfun(mfun(preds_all));
accuracy_son_psth_sem = cfun(std(preds_all,[],1,'omitnan')/sqrt(nfold));
%tests
d_sh = d_son_shuffle;
d_ofc = d_son_ofc;
accuracy_pred_means_ofc = mean(d_ofc.preds_all_ofc,2,'omitnan');
accuracy_pred_means_shuffle = mean(d_sh.preds_all_shuffle,2,'omitnan');
[pgood_mlb_son, pv_son_shuff, pv_son_ofc] = sigtests(preds_all,accuracy_pred_means_ofc,accuracy_pred_means_shuffle,output.xvec, alpha);


d = d_son2;
preds_all = zeros(nfold,ntx);
for m = 1:nfold
    for j = 1:ntx
        preds_all(m,j) = sum(d.preds_CV(m,j,:)==d.true_CV(m,j,:))/numel(d.true_CV(m,j,:));
    end
end
accuracy_son_psth_mean2 = cfun(mfun(preds_all));
accuracy_son_psth_sem2 = cfun(std(preds_all,[],1,'omitnan')/sqrt(nfold));
%tests
d_sh = d_son_shuffle2;
d_ofc = d_son_ofc2;
accuracy_pred_means_ofc = mean(d_ofc.preds_all_ofc,2,'omitnan');
accuracy_pred_means_shuffle = mean(d_sh.preds_all_shuffle,2,'omitnan');
[pgood_mlb2_son, pv_son_shuff2, pv_son_ofc2] = sigtests(preds_all,accuracy_pred_means_ofc,accuracy_pred_means_shuffle,output.xvec, alpha);




%%  plot average and sem predictive performance over time. both decoders single plot

%fac = 1.96*sqrt(nfold); %is using 95% CI for ofc, use this
fac = 1; % if using confidence sem, use this

%y lims
%yl = [0,1];
yl = [0.25,0.5];
xvec = output.xvec;

figure(255)
clf

% COFF
subplot(2,3,1)
hold on
shadedErrorBar(xvec_ofc,accuracy_coff_psth_mean, accuracy_coff_psth_sem,'lineprops',{'color','red'})
%shadedErrorBar(xvec_ofc,d_coff_ofc.accuracy_ofc_mean, fac*d_coff_ofc.accuracy_ofc_sem,'lineprops',{'color','blue'})
shadedErrorBar(xvec_ofc,cfun(d_coff_ofc.accuracy_ofc_median), cfun(fac*d_coff_ofc.accuracy_ofc_sem),'lineprops',{'color','blue'})
shadedErrorBar(output.xvec,cfun(d_coff_shuffle.accuracy_shuffle_mean), cfun(1.96*d_coff_shuffle.accuracy_shuffle_std),'lineprops',{'color','black'})
plot(xvec(find(pv_coff_shuff < alpha)),0.29*ones(size(find(pv_coff_shuff < alpha))),'*', 'linewidth', 0.5, 'color','k')
% plot(xvec(find(pv_coff_ofc < alpha)),0.27*ones(size(find(pv_coff_ofc < alpha))),'*', 'linewidth', 0.5, 'color','b')
title('MLB, epoch = COFF')
xlabel('Time from tone (s)')
ylabel('Accuracy (correct)')
legend('OFC-CPi', 'OFC','shuffle')
ylim(yl)
xline(0,'--k')
%xlim([-4,4])
xlim([0,3])
set(gca,'fontsize',15)


subplot(2,3,4)
hold on
shadedErrorBar(xvec_ofc,accuracy_coff_psth_mean2, accuracy_coff_psth_sem2,'lineprops',{'color','red'})
%shadedErrorBar(xvec_ofc,d_coff_ofc2.accuracy_ofc_mean, fac*d_coff_ofc2.accuracy_ofc_sem,'lineprops',{'color','blue'})
shadedErrorBar(xvec_ofc,cfun(d_coff_ofc2.accuracy_ofc_median), cfun(fac*d_coff_ofc2.accuracy_ofc_sem),'lineprops',{'color','blue'})
shadedErrorBar(output.xvec,cfun(d_coff_shuffle2.accuracy_shuffle_mean), cfun(1.96*d_coff_shuffle2.accuracy_shuffle_std),'lineprops',{'color','black'})
plot(xvec(find(pv_coff_shuff2 < alpha)),0.29*ones(size(find(pv_coff_shuff2 < alpha))),'*', 'linewidth', 0.5, 'color','k')
% plot(xvec(find(pv_coff_ofc2 < alpha)),0.27*ones(size(find(pv_coff_ofc2 < alpha))),'*', 'linewidth', 0.5, 'color','b')
title('MLB2, epoch = COFF')
xlabel('Time from tone (s)')
ylabel('accuracy')
ylim(yl)
xline(0,'--k')
%xlim([-4,4])
xlim([0,3])
set(gca,'fontsize',15)


% SON: incomplete
subplot(2,3,2)
hold on
shadedErrorBar(xvec_ofc,accuracy_son_psth_mean, accuracy_son_psth_sem,'lineprops',{'color','red'})
%shadedErrorBar(xvec_ofc,d_son_ofc.accuracy_ofc_mean, fac*d_son_ofc.accuracy_ofc_sem,'lineprops',{'color','blue'})
shadedErrorBar(xvec_ofc,cfun(d_son_ofc.accuracy_ofc_median), cfun(fac*d_son_ofc.accuracy_ofc_sem),'lineprops',{'color','blue'})
shadedErrorBar(output.xvec,cfun(d_son_shuffle.accuracy_shuffle_mean), cfun(1.96*d_son_shuffle.accuracy_shuffle_std),'lineprops',{'color','black'})
plot(xvec(find(pv_son_shuff < alpha)),0.29*ones(size(find(pv_son_shuff < alpha))),'*', 'linewidth', 0.5, 'color','k')
% plot(xvec(find(pv_son_ofc < alpha)),0.27*ones(size(find(pv_son_ofc < alpha))),'*', 'linewidth', 0.5, 'color','b')
title('MLB, epoch = SON')
xlabel('Time from side on (s)')
ylabel('Accuracy (correct)')
ylim(yl)
xline(0,'--k')
xlim([-4,4])
xlim([-1,3])
set(gca,'fontsize',15)


subplot(2,3,5)
hold on
shadedErrorBar(xvec_ofc,accuracy_son_psth_mean2, accuracy_son_psth_sem2,'lineprops',{'color','red'})
%shadedErrorBar(xvec_ofc,d_son_ofc2.accuracy_ofc_mean, fac*d_son_ofc2.accuracy_ofc_sem,'lineprops',{'color','blue'})
shadedErrorBar(xvec_ofc,cfun(d_son_ofc2.accuracy_ofc_median), cfun(fac*d_son_ofc2.accuracy_ofc_sem),'lineprops',{'color','blue'})
shadedErrorBar(output.xvec,cfun(d_son_shuffle2.accuracy_shuffle_mean), cfun(1.96*d_son_shuffle2.accuracy_shuffle_std),'lineprops',{'color','black'})
plot(xvec(find(pv_son_shuff2 < alpha)),0.29*ones(size(find(pv_son_shuff2 < alpha))),'*', 'linewidth', 0.5, 'color','k')
% plot(xvec(find(pv_son_ofc2 < alpha)),0.27*ones(size(find(pv_son_ofc2 < alpha))),'*', 'linewidth', 0.5, 'color','b')
title('MLB2, epoch = SON')
xlabel('Time from side on (s)')
ylabel('accuracy')
ylim(yl)
xlim([-4,4])
xlim([-1,3])
xline(0,'--k')
set(gca,'fontsize',15)


% reward
subplot(2,3,3)
hold on
shadedErrorBar(xvec_ofc,accuracy_reward_psth_mean, accuracy_reward_psth_sem,'lineprops',{'color','red'})
%shadedErrorBar(xvec_ofc,d_reward_ofc.accuracy_ofc_mean, fac*d_reward_ofc.accuracy_ofc_sem,'lineprops',{'color','blue'})
shadedErrorBar(xvec_ofc,cfun(d_reward_ofc.accuracy_ofc_median), cfun(fac*d_reward_ofc.accuracy_ofc_sem),'lineprops',{'color','blue'})
shadedErrorBar(output.xvec,cfun(d_reward_shuffle.accuracy_shuffle_mean), cfun(1.96*d_reward_shuffle.accuracy_shuffle_std),'lineprops',{'color','black'})
plot(xvec(find(pv_reward_shuff < alpha)),0.29*ones(size(find(pv_reward_shuff < alpha))),'*', 'linewidth', 0.5, 'color','k')
% plot(xvec(find(pv_reward_ofc < alpha)),0.27*ones(size(find(pv_reward_ofc < alpha))),'*', 'linewidth', 0.5, 'color','b')
title('MLB, epoch = reward')
xlabel('Time from reward (s)')
ylabel('accuracy')
ylim(yl)
xline(0,'--k')
xlim([-4,4])
xlim([-1,3])
set(gca,'fontsize',15)


subplot(2,3,6)
hold on
shadedErrorBar(xvec_ofc,accuracy_reward_psth_mean2, accuracy_reward_psth_sem2,'lineprops',{'color','red'})
%shadedErrorBar(xvec_ofc,d_reward_ofc2.accuracy_ofc_mean, fac*d_reward_ofc2.accuracy_ofc_sem,'lineprops',{'color','blue'})
shadedErrorBar(xvec_ofc,cfun(d_reward_ofc2.accuracy_ofc_median), cfun(fac*d_reward_ofc2.accuracy_ofc_sem),'lineprops',{'color','blue'})
shadedErrorBar(output.xvec,cfun(d_reward_shuffle2.accuracy_shuffle_mean), cfun(1.96*d_reward_shuffle2.accuracy_shuffle_std),'lineprops',{'color','black'})
plot(xvec(find(pv_reward_shuff2 < alpha)),0.29*ones(size(find(pv_reward_shuff2 < alpha))),'*', 'linewidth', 0.5, 'color','k')
% plot(xvec(find(pv_reward_ofc2 < alpha)),0.27*ones(size(find(pv_reward_ofc2 < alpha))),'*', 'linewidth', 0.5, 'color','b')
title('MLB2, epoch = reward')
xlabel('Time reward (s)')
ylabel('accuracy')
ylim(yl)
xline(0,'--k')
xlim([-4,4])
xlim([-1,3])
set(gca,'fontsize',15)


sgtitle(strcat('subtraction type = ',decodertype,', is stratified: ',string(is_stratified)),'fontsize',20)

for s = 1:6
    subplot(2,3,s)
    xlim([-1 2])
    set(gca,'TickDir','out','Box','off')
end
subplot(2,3,1)
xlim([0 2])
subplot(2,3,4)
xlim([0 2])


%% aggregate and save pvalues to csv


tvec = xvec;

D = [tvec; ...
         pv_coff_shuff; pv_coff_ofc; ...
         pv_coff_shuff2; pv_coff_ofc2;...
         pv_son_shuff; pv_son_ofc; ...
         pv_son_shuff2; pv_son_ofc2;...
         pv_reward_shuff; pv_reward_ofc; ...
         pv_reward_shuff2; pv_reward_ofc2]';

writematrix(D, strcat(savedir,'pvalues.csv'))

%% dignificance test code

function [pgood, pvals_shuff, pvals_ofc] = sigtests(dat,dat_ofc,dat_shuffle,xvec, alpha)
% runs significance testing for 
    xvec_ofc = -4:0.2:4;
    

    nt = numel(xvec);
    pvals_ofc = nan(1,nt);
    pvals_shuff = nan(1,nt);
    for j = 1:numel(xvec)
    
        [~,idx_ofc] = min(abs(xvec_ofc-xvec(j)));
        %[~,idx] = min(abs(xvec-xvec(j)));
        idx = j;
        dat_ofc_j = dat_ofc(:,idx_ofc);
        dat_shuffle_j = dat_shuffle(:,idx);
        dat_j = dat(:,idx);
        
        %pvals_ofc(j) = signrank(dat_ofc_j,dat_j);
        %%[~,pvals_shuff(j)] = ttest(dat_shuffle_j,dat_j);
        %pvals_shuff(j) = signrank(dat_shuffle_j,dat_j);

        pvals_ofc(j) = ranksum(dat_ofc_j,dat_j);
        %[~,pvals_shuff(j)] = ttest(dat_shuffle_j,dat_j);
        pvals_shuff(j) = ranksum(dat_shuffle_j,dat_j);

    end
    
    %pvals_ofc(pvals_ofc > alpha) = nan;
    %pvals_shuff(pvals_shuff > alpha) = nan;
    %pgood = find(~isnan(pvals_shuff & pvals_ofc));
    pgood = find(pvals_shuff < alpha & pvals_ofc < alpha);

end


