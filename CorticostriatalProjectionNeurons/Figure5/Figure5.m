function Figure5(datapath,codepath, savedir)
% datapath = 'C:\Users\mld9131\Documents\GitHub\Papers\MLD\Figure5\Data';
% codepath = 'C:\Users\mld9131\Documents\GitHub\Papers\MLD\Figure5\Code'
% savepath =  Z:\david\maggie\data\

%% A. Schematic of the simplex
% add path and load DS projection data

%location of analysis codebase
% codepath = '~/david/maggie';

% add code paths
% addpath(genpath(codepath))
% addpath(genpath(codepath + "david"))
% addpath(genpath(codepath + "maggie"))

%for saving pds correctly so illustrator can modify them
set(0, 'DefaultFigureRenderer', 'painters');

% DS projection neuron data location
%datapath = "/Users/dhocker/projects/dynamics/data/maggie/";
% datpath = "Z:\david\maggie\Most_Likely_block_analysis\";
fname = fullfile(datapath, "DS_projection_neurons_non-Stimulated.mat");
b = load(fname);

%% run bayes model for each session and get MLB, MLB2
% for every cell (only 60 now, and 23 sessions, so not prohibitive), do bayes calc
% bayes params for all sessions
kappa_mi = 1; %.5;  
kappa_hi = 1.5; %1.2; %0.7; 
kappa_lo = 0.5;%0.8; %.3; 
D = 0.5; 
lambda =1;
noise = 80;
params = [kappa_mi, kappa_hi, kappa_lo, D, lambda];
nsess  = numel(b.all_S_DS);

Beliefs_all = [];
wait_time_all = [];
Beliefs_bysess = {};
for j = 1:nsess
    disp(j)
    %get session data
    S_j = b.all_S_DS{j};

    % this struct is almost the correct format for an A struct, but needs
    % a little renaming to run. 
    S_j.ntrials = [numel(S_j.Block)]; % only 1 session in this struct
    S_j.reward = S_j.RewardAmount;
    S_j.prob_catch = S_j.ProbCatch;

    [~, wait_time, ~, Belief, ~, ~] = GenerateSynthData_Bayes(params, S_j,'logn', true,noise);
    Beliefs_bysess{j} = Belief;
    Beliefs_all = [Beliefs_all, Belief];
    wait_time_all = [wait_time_all, wait_time'];

end

%remove nans
Beliefs_m = Beliefs_all(1,:);
Beliefs_h = Beliefs_all(2,:);
Beliefs_l = Beliefs_all(3,:);
Beliefs_m(isnan(wait_time_all)) = [];
Beliefs_h(isnan(wait_time_all)) = [];
Beliefs_l(isnan(wait_time_all)) = [];

Beliefs_scrubbed = [Beliefs_m;Beliefs_h; Beliefs_l];
    
%calculate mostlikely block (MLB) and 2nd most likely block (MLB2)
[B,I] = sort(Beliefs_scrubbed,1,'descend');
mlb = I(1,:);
mlb2 = I(2,:);

%% plot simplex
figure('Color','w')
hold on

grid
xlabel('p mix')
ylabel('p high')
zlabel('p low')
title('belief distribution')
set(gca,'fontsize',15)

 
x = [ 1, 0.5, 1/3, 0.5];
y = [ 0, 0.5, 1/3, 0];
z = [ 0, 0,   1/3, 0.5];
fill3(x,y,z,'k','facealpha',0.2)

%high
x = [ 0, 0.5, 1/3, 0.0];
y = [ 1, 0.5, 1/3, 0.5];
z = [ 0, 0,   1/3, 0.5];
fill3(x,y,z,'r','facealpha',0.2)

%low
x = [ 0, 0.0, 1/3, 0.5];
y = [ 0, 0.5, 1/3, 0.0];
z = [ 1, 0.5,   1/3, 0.5];
fill3(x,y,z,'b','facealpha',0.2)
view(69,25)

% plot beliefs, color coded by MLB2. give noise jitter to location for vis.
ntrial = size(Beliefs_scrubbed,2);
noise = 0.005*randn(3,ntrial);
scatter3(Beliefs_scrubbed(1,mlb2==1)+noise(1,mlb2==1), Beliefs_scrubbed(2,mlb2==1)+noise(2,mlb2==1), Beliefs_scrubbed(3,mlb2==1)+noise(3,mlb2==1),30,'o','markerfacecolor','k','markeredgecolor','k','markerfacealpha',0.2,'markeredgealpha',0.2)
scatter3(Beliefs_scrubbed(1,mlb2==2)+noise(1,mlb2==2), Beliefs_scrubbed(2,mlb2==2)+noise(2,mlb2==2), Beliefs_scrubbed(3,mlb2==2)+noise(3,mlb2==2),30,'o','markerfacecolor','r','markeredgecolor','r','markerfacealpha',0.2,'markeredgealpha',0.2)
scatter3(Beliefs_scrubbed(1,mlb2==3)+noise(1,mlb2==3), Beliefs_scrubbed(2,mlb2==3)+noise(2,mlb2==3), Beliefs_scrubbed(3,mlb2==3)+noise(3,mlb2==3),30,'o','markerfacecolor','b','markeredgecolor','b','markerfacealpha',0.2,'markeredgealpha',0.2)

ylim([-0.0,1])
xlim([-0.0,1])
zlim([-0.0,1])


%% B. 20uL encoding when p(mix)>0.8 and p(high) or p(low) <0.2
datapath2 = 'C:\Users\mld9131\Documents\GitHub\Papers\MLD\Figure2\Data';
load(fullfile(datapath2,'OFC-CPi_projection_neurons_non-Stimulated.mat'));
% load(fullfile('\\constantinoplelab.cns.nyu.edu\server2\PhysiologyData\Maggie\Chronic_implant\Npxl\Optotagged_cells','OFC-DS_projection_neurons_non-Stimulated.mat'))
[pval,perlow] = uncertainty_exploratory(all_SU_DS, all_index_DS, all_S_DS);

times = -1:0.05:3;
xvec = all_SU_DS{1, 1}.xvec.CON;
t = xvec>=-1 & xvec<=3;
pval2 = pval(t);
T = table(times, pval2,'VariableNames',{'Time','p value'});
% writetable(T,'Z:\Maggie\Papers\Physiology\Extended_data_tables\Fig5d_20uL-2MLB.xlsx')

%% C & D:  Decoding Most likely block and Second most likely block
% 'Z:\david\maggie\code\study_20260609_figgen_allepoch.m'
%% load code path and data.

%location of analysis codebase
% codepath = '~/projects/constantinoplelab/Analysis/';
% 
% % add code paths
% addpath(genpath(codepath))
% addpath(genpath(codepath + "david"))
% addpath(genpath(codepath + "maggie"))

%for saving pds correctly so illustrator can modify them
set(0, 'DefaultFigureRenderer', 'painters');

% DS projection neuron data location
% datapath = '\\constantinoplelab.cns.nyu.edu\server\david\maggie\data\'; % "/Users/dhocker/projects/dynamics/data/maggie/";
% savedir = 'Z:\david\maggie\data\'; %"/Users/dhocker/projects/dynamics/results/maggie/";

fname = fullfile(datapath,"DS_projection_neurons_non-Stimulated.mat");
b = load(fname);

% parse the DS projections to get beliefs
usetorch = false;
% output = dSTR_decode_parsedata(usetorch);
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

%% calculate the accuracy for the projection neurons

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
[pgood_mlb_reward, pv_reward_shuff, pv_reward_ofc] = sigtests(preds_all,accuracy_pred_means_ofc,accuracy_pred_means_shuffle,xvec_ofc, alpha);



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

% figure
nexttile(2)
% COFF
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

xlim([0,2])
set(gca,'fontsize',15,'TickDir','out','Box','off')


nexttile(6)
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

xlim([0,2])
set(gca,'fontsize',15,'TickDir','out','Box','off')


% SON
nexttile(3)
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
   xlim([-1 2])
set(gca,'fontsize',15,'TickDir','out','Box','off')


nexttile(7)
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
   xlim([-1 2])
xline(0,'--k')
set(gca,'fontsize',15,'TickDir','out','Box','off')


% reward
nexttile(4)

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
   xlim([-1 2])
set(gca,'fontsize',15,'TickDir','out','Box','off')


nexttile(8)
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
xlim([-1 2])
set(gca,'fontsize',15,'TickDir','out','Box','off')



sgtitle(strcat('subtraction type = ',decodertype,', is stratified: ',string(is_stratified)),'fontsize',20)


%% E: Schematic of Preferred encoding
%% F: Preferred v non--preferred block encoding

% using just OFC-CPI data this will look for the preferred block at COff.
% load the data
datpath = 'C:\Users\mld9131\Documents\GitHub\Papers\MLD\Figure2\Data';
load(fullfile(datpath,'OFC-CPi_projection_neurons_non-Stimulated.mat'));

% initialize
xvec_200 = -4:0.2:8;
xvec_100 = -4:0.1:8;
xvec = -4:0.05:8;
vcolors = [.6 .2 .9; .1 .5 .9; 0 .4 0; 1 .65 0; 1 0 0]; %need to play with colors to be color blind friendly
alignto = {'COFF', 'SON','SOFF','Rew','Opt','CON'};

%generate z-score fr
nanzscore = @(x)(x - mean(x,2, 'omitnan'))./std(x,0,2, 'omitnan');
for k = 1:4
for i = 1:length(DLS)
    all_SU_z{1,DLS(i)}.hmat.(string(alignto(k))) = nanzscore(all_SU_DS{1,DLS(i)}.hmat.(string(alignto(k))));
end
end

for i = 1:length(DLS)
    S = all_S_DS{all_index_DS{i,4}};
    TL = S.Block==3; % & S.RewardAmount<20;
    TH = S.Block==2; % & S.RewardAmount>20;
    TM = S.Block==1; %

    win = [0 1];
    T = find(xvec>=win(1),1):find(xvec<=win(2),1,'last');

   
    % Determine whether high or low blocks are the preferred block
    FRL = mean(mean(all_SU_DS{i}.hmat.COFF(TL,T),'omitnan'));
    FRH =  mean(mean(all_SU_DS{i}.hmat.COFF(TH,T),'omitnan'));
    FRM = mean(mean(all_SU_DS{i}.hmat.COFF(TM,T),'omitnan'));
    reward = convertreward(S.RewardAmount);
    for k = 1:4
        for rew = 1:5
            if rew<3
                if FRL>FRM 
                    preftrials = reward==rew & S.Block==3;
                    nonpreftrials = reward==rew & S.Block==1;
                else
                    preftrials = reward==rew & S.Block==1;
                    nonpreftrials = reward==rew & S.Block==3;
                end

            elseif rew>3
                if FRM>FRH
                    preftrials = reward==rew & S.Block==1;
                    nonpreftrials = reward==rew & S.Block==2;
                else
                    preftrials = reward==rew & S.Block==2;
                    nonpreftrials = reward==rew & S.Block==1;
                end

            else
                if FRL>FRH
                    preftrials = reward==rew & S.Block==3;
                    nonpreftrials = reward==rew & S.Block==2;
                else
                    preftrials = reward==rew & S.Block==2;
                    nonpreftrials = reward==rew & S.Block==3;
                end
                mixtrials = reward==rew & S.Block==1;
                FRMix.(alignto{k}){rew}(i,:) = mean(all_SU_DS{i}.hmat.(alignto{k})(mixtrials,:),'omitnan');

            end

            FRpref.(alignto{k}){rew}(i,:) = mean(all_SU_DS{i}.hmat.(alignto{k})(preftrials,:),'omitnan');
            FRnonpref.(alignto{k}){rew}(i,:) = mean(all_SU_DS{i}.hmat.(alignto{k})(nonpreftrials,:),'omitnan');

            FRprefz.(alignto{k}){rew}(i,:) = mean(all_SU_z{i}.hmat.(alignto{k})(preftrials,:),'omitnan');
            FRnonprefz.(alignto{k}){rew}(i,:) = mean(all_SU_z{i}.hmat.(alignto{k})(nonpreftrials,:),'omitnan');
        end
    end
end


%% Plot the response
z = 0;
win =[-1 2.5];
T = find(xvec>=win(1),1):find(xvec<=win(2),1,'last');
   for k = [1,2,4]
        nexttile
        if z==1
            % shadedErrorBar(xvec(T),mean(FRprefz.(alignto{k}){rew}(:,T),'omitnan'),sem(FRprefz.(alignto{k}){rew}(:,T),'omitnan'),'Lineprops',{'Color','k'});
            plot(xvec(T),mean(FRprefz.(alignto{k}){rew}(:,T),'omitnan'),'Color',[0 0 0]);
            hold on
            % shadedErrorBar(xvec(T),mean(FRnonpref.(alignto{k}){rew}(:,T),'omitnan'),sem(FRnonprefz.(alignto{k}){rew}(:,T),'omitnan'),'Lineprops',{'Color',[0.35 0.35 .35]});
            plot(xvec(T),mean(FRnonprefz.(alignto{k}){rew}(:,T),'omitnan'),'Color', [0.5 0.5 0.5]); %[.9 0.7 .9]);

        else
            % shadedErrorBar(xvec(T),mean(FRpref.(alignto{k}){rew}(:,T),'omitnan'),sem(FRpref.(alignto{k}){rew}(:,T),'omitnan'),'Lineprops',{'Color','k'});
            plot(xvec(T),mean(FRpref.(alignto{k}){rew}(:,T),'omitnan'),'Color',[0 0 0]);
            hold on
            % shadedErrorBar(xvec(T),mean(FRnonpref.(alignto{k}){rew}(:,T),'omitnan'),sem(FRnonpref.(alignto{k}){rew}(:,T),'omitnan'),'Lineprops',{'Color',[0.35 0.35 .35]});
            plot(xvec(T),mean(FRnonpref.(alignto{k}){rew}(:,T),'omitnan'),'Color', [0.5 0.5 0.5]); %[.9 0.7 .9]);

            % % shadedErrorBar(xvec(T),mean(FRMix.(alignto{k}){rew}(:,T),'omitnan'),sem(FRMix.(alignto{k}){rew}(:,T),'omitnan'),'Lineprops',{'Color',[0.0 0.5 .6]});
            % plot(xvec(T),mean(FRMix.(alignto{k}){rew}(:,T),'omitnan'),'Color',[.7 0.3 .7]);

            ylim([3.5 10])
        end

        set(gca,'Box','off','TickDir','out')
        xlabel('Time from event (s)');xline(0,'--');xlim(win);
        title(alignto{k})
   end
legend({'Preferred','NonPreferred'},'location','southeast')
ylabel('Firing rate (Hz)')
vol = [5,10,20,40,80];
sgtitle(strcat("Volume = ",string(vol(rew)),'uL'))

end



%% significance test code

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


