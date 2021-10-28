%demo for running the clustering pipeline and the gap statistic

%% create the intial building blocks for features
ops = struct();

%set up a form often used: 90% bootstrap sampling, with repalcement.
%z-scored
ops.nsamps = 100; %how many bootstrap for ARI samples will be made?
ops.usezscore = true; %z-scored responses as feature space
ops.bootstrap_proportion = 0.9; %proportion of samples in bootstrap
ops.bootstrap_replacement = true; % with replacement= = true
ops.normfet = 0.9; %is z-score is false, should features be normalzied by max response?
ops.savedir = '/Users/dhocker/projects/ofc/results/demo/'; %where will data be stored
ops.datadir = '~/projects/ofc/data/published/';

[fname_base] = sourceFeatures_ofc(ops); %source features from basic aspects of data

ops.fname_base = fname_base; %the save name for the source features


%% or just load a previous set
%ops = struct();
%ops.savedir = '/Users/dhocker/projects/ofc/results/demo/'; %where will data be stored
%ops.datadir = '~/projects/ofc/data/published/';

%fname_base = 'ofcFeaturebase_Z_bs90Rep_ns_100.mat';
%ops.psthfile = 'psth_conditional_all.mat';
%ops.fname_basemodel = 'ofcFeaturebase_models.mat';


%f = load(strcat(ops.datadir,fname_base));
%p = load(strcat(ops.datadir,ops.psthfile));
%fmod = load(strcat(ops.datadir,ops.fname_basemodel));


%% choose a specific feature type and build that feature space

disp('creating specific feature space representation')

%PSTH clustering

ops.fname_basemodel = 'ofcFeaturebase_models.mat';
ops.featureClass = 'PSTH'; %feature space class, def. in createFeatures_ofc.m
ops.featureType = 1; %specific feature space type in featureClass
ops.preprocType = 'pca'; %what preprocessing to use?
ops.preprocCutoff = 0.95; %95% variance cutoff 


%conditional clustering
%{
ops.fname_basemodel = 'ofcFeaturebase_models.mat';
ops.featureClass = 'singleConditional';
ops.featureType = 12;
ops.fname_base = fname_base;
ops.preprocType = 'pca';
ops.preprocCutoff = 0.95; %95% variance cutoff 
%}


[fet_proc,fet,fet_sem,fet_legend,xvec] =...
    createFeatures_ofc(ops);

ops.xvec = xvec;

%write to file
savename = fullfile(ops.savedir,'features.mat');
save(savename,'fet_proc','fet','fet_sem','fet_legend','ops');

%% calculte PAIRS

fet_proc = f.fet_proc(:,:,1);
kval = 3; %for cond clustering. k=8. for psth clustering: k=3
nset = 100; %for real comp. use 10k
saveops = struct();
saveops.saveref = false; %save the reference data?
saveops.loadref = true; %load the reference data?
%saveops.savename = strcat(ops.savedir,'pairs_ref.mat');
saveops.loadname = strcat(ops.savedir,'pairs_ref.mat');

[pairs,p_pairs,pairs_ref,meanvals_dat,meanvals] = ...
        PAIRS(fet_proc,kval,nset,saveops);
    
    
%plot results
graphgen_PAIRS(pairs,pairs_ref,meanvals_dat,meanvals);


%% call the gap statistic, then plot

f = load(savename);
fet = f.fet(:,:,1);
fet_proc = f.fet_proc(:,:,1); %preprocessed one

%nn = size(fet,1);
clustvec = 2:15;


nsamp = 100; %number of null dist samples %5k to 10k is best. small here for speed

% calculate gap score
disp('running gap statistic')
eva = evalclusters(fet,'kmeans','gap','KList',clustvec,'SearchMethod','firstMaxSE','B',nsamp);
labels = eva.OptimalY; %cluster labels
N = eva.OptimalK; %number of located clusters

figure(93)
clf
plot(eva)
vline(N,'k')
title('gap statistic results (error bars = 1SEM)')
set(gca,'fontsize',15)

%save
savename_gap = fullfile([ops.savedir,'gap.mat']);
save(strcat(savename_gap,'eva_kmeans.mat'),'eva','labels','N')






