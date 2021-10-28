%fits a sample GLM 
%written by dlh, Oct. 2021

%runtime is about 15 mintes on a laptop, with hyperparameter search incl.

%% add paths and set up load and save directories

codedir = '/Users/dhocker/projects/ofc';
savedir = fullfile('/Users/dhocker/projects/ofc/results/'); 
datadir = '/Users/dhocker/projects/ofc/data/published/units/';

addpath(genpath(codedir));


%% set up glm parameters for a sample file

namebase = 'fit_stim2020903_N'; %name prefix for saved glm file
popname = 'pop_stim2020903'; %the savedir suffix %TODO: needed?
file_id = 22; %TODO: find best neuron

dataname = strcat([datadir,'unit_',num2str(file_id),'.mat']);


%parameters proper. see glm_ofcOps() for all of the settings
ops =  glm_ofcOps(); %common ofc parameters. easy loading script    
ops.fname = strcat([savedir,namebase,num2str(file_id)]); %save name
ops.dataname = dataname; %original data file

%how will optimization be handled
ops.useHyperParamGrid = true; %do grid search for hyperparameters?
ops.kfold = 5; %how many kfolds for cross val? includes 1 for test set

%decide on stimulus inputs. make sure ops.nkern is set correctly, and 
%corresponds to number of covariates in ops.stimtype, defined in loadGlmDat_ofc.
ops.stimtype = '20200903'; %set the types of inputs in the model
ops.nkern = 20; %how many kernels in model from ops.stimtype.

%set the dimensionality and support for each of the nkern kernels to match
%nkern
ls = ceil(4/ops.dt); %length of kernel, in timesteps
ops.ls = ls*ones(1,ops.nkern); %
ops.dimst = 9*ones(1,ops.nkern); %how many basis functions per kernel
ops.dimbcst = 1*ones(1,ops.nkern); %how many boxcar functions in basis

%% fit the glm

fun_glmfit(ops);

%% visualize kernels, comparison to data

%load model
f = load(ops.fname);

basefuns = f.stim.basefuns{2}; %the basis functions
stimnames = f.stim.stimlegend; %names of the kernels
tvec = f.ops.dt:f.ops.dt:f.ops.ls(1)*f.ops.dt; %all ls are the same for this model


%which k fold was best based on validation set
[~,kind] = min(f.NLL.val);

params = f.wML{kind};
params_cov = f.wML_cov{kind};

nkern = numel(basefuns);
[dim,ntkern] = size(basefuns{1});
kernels = cell(nkern,1);

figure(20);
ind = 0;
for k = 1:nkern
    subplot(4,ceil(nkern/4),k)
    kernels{k} = params(ind+1:ind+dim)*basefuns{k};  
    plot(tvec,kernels{k},'linewidth',2);
    title(stimnames{k})
    
    ind = ind + dim;
end


%% compare model to data


%some parameters to simulate and compare

trialmask = f.ops.test_ids; %only compare to test data
%trialmask = true(size(f.ops.test_ids)); %compare to all data

alignment = 'start';
tmin = -2; %start stimulation time from alignment event
tmax = 4; %finish simulation time from alignment event



%parse up trials and simulate from model
[rate_dat,lambda,trial_idx_clip,lambda_var] = cleantrials(f,trialmask,alignment,tmin,tmax,ops);


condmask = true(size(trial_idx_clip));

tvec = tmin:ops.dt:tmax;

psth_dat = nanmean(rate_dat(condmask,:));
psth_dat_sem = nanstd(rate_dat(condmask,:))/sqrt(sum(condmask));
psth_mod = nanmean(lambda(condmask,:));
psth_mod_sem =  nanstd(lambda(condmask,:))/sqrt(sum(condmask));



figure(147)
clf
shadedErrorBar(tvec,psth_dat,psth_dat_sem,...
    'lineprops',{'color','k','linewidth',2});
hold on
shadedErrorBar(tvec,psth_mod,psth_mod_sem,...
    'lineprops',{'color','r','linewidth',2});
xlabel(strcat('time to ',alignment,' (s)'))
ylabel('rate (Hz)')
title(strcat('model comparison, Neuron ',num2str(file_id)))
vline(0,'k')
legend('data','model')
set(gca,'fontsize',15)



