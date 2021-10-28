%calcualtes CPD for an example cell

%% choose a neuron and its glm fit, as well as place to save

file_idx = 22;

datadir = '~/projects/ofc/data/published/units/'; %where raw data lives
datadir_agg = '~/projects/ofc/data/published/'; %where aggregated data lives
savedir = '~/projects/ofc/results/demo/'; %where CPD should be saved

fname_glm = strcat(savedir,'fit_stim2020903_N',num2str(file_idx),'.mat');


%% load the A struct of concatenated data 
%to help with creating trial contingency masks

D = load(strcat(datadir_agg,'concatdata_ofc_noOptOut.mat'));
A = D.A;
clear D

%% set up helper cells to organize how covariates will be omitted/removed
% for the reduced model in CPD

%sets of kernels to be be averaged or fully omitted for a given CPD calc.
groups_cpd = {{'prevWin','prevLoss','prevOptOut','prevWin (ac)','prevLoss (ac)','prevOptOut (ac)'},...
            {'Left','Right'},...
            {'win','loss'},...
            {'risky (ac)','safe (ac)'},...
            {'LFlash','RFlash'},...
            {'LClick','RClick'}};

%what event should be used to align each CPD calc
alignlist_cpd =  {'start','leavecpoke','choice','choice','start','start'};
        
%naming convention for each CPD calculation. used for file saving
groupnames_cpd = {'rewardHistory','LR','reward','riskySafe',...
    'flash','click'};  

%conditional masks for each type. used to weight each sample based on trial type    
condgroups = {{'prevWin','prevLoss','prevOptOut'},...
            {'left','right'},...
            {'win','loss'},...
            {'risky','safe'},...
            '',''}; %clicks and flashes use all trials, so no cond. mask 
        
%generate trial masks for each type of trial contengincy
condmask = cell(numel(groupnames_cpd),1);
for j = 1:numel(groupnames_cpd)    
    %the common case, sets of binary conditions, one per trial, and easily
    %distinguished as a case from condmaskLibrary with a string handle
    if iscell(condgroups{j})
        cj = cell(numel(condgroups{j}),1);
        for k = 1:numel(condgroups{j})
            cj{k} = condmaskLibrary(condgroups{j}{k} ,file_idx,A,struct('removeOptOut',true));
        end
    elseif strcmp(condgroups{j},'') % a null case. no weighting needed. set to all trials
        cj = condmaskLibrary('prevWin' ,file_idx,A,struct('removeOptOut',true)); %chose a random condition
        cj = true(size(cj)); %set all trials true
    else %case for scalar condition with single 'string' handle. e.g., 'vol'
        cj = condmaskLibrary(condgroups{j} ,file_idx,A,struct('removeOptOut',true));
    end
    
    condmask{j} = cj;
end

%time ranges for each alignment type, in seconds from alignment event
tminvec = [-2,-4,-4];
tmaxvec = [4,4,4];
alcell = {'start','leavecpoke','choice'};


%% set the options for the code

j = 1; %index in the helper lists above. which cpd calc to perform


disp(strcat(['CPD TYPE: ',groupnames_cpd{j}]))
ops = struct();
ops.datadir = datadir; %where raw data lives. TODO delete?
ops.fname_glm = fname_glm; %filename of glm  fit being used for calc.
ops.alignment = alignlist_cpd{j}; %what even is used for alignment 
ops.masktype = 'all'; %'test' data or 'all' data? 
ops.tmin = tminvec(strcmp(alignlist_cpd{j},alcell)); %trial start time (s)
ops.tmax = tmaxvec(strcmp(alignlist_cpd{j},alcell)); %trial end time (s)
ops.simtype = 'glm'; %how to simulate form GLM
ops.cpdgroup = groups_cpd{j};  %which kernels to average/omit
ops.cpdname = groupnames_cpd{j}; %name of CPD being calculated. TODO remove
ops.shufftype = 'behavior'; %'behavior'= averaging kernels. 'omit' = remove
ops.condmask = condmask{j};%contingency masks for trial balancing

%%
[Vt,Vtshuff,shuffmean,shuffstd] = CPD_noOpt(file_idx,ops);

%calculate the significance region:
%get shuffle-mean substracted, significant CPD only
[h,x] = histcounts(Vtshuff,20,'Normalization','cdf');
g95 = find(h > 0.975);
shuffsig = x(g95(1)); %first instance
Vtsig = Vt-shuffmean;
Vtsig(Vtsig < shuffsig) = 0;
    
%% save

savename = strcat(savedir,groupnames_cpd{j},'_',num2str(file_idx),'.mat');

save(savename,'Vtsig','Vt','Vtshuff','shuffmean','shuffstd','shuffsig')
    

%% plot

figure(110)
clf
tmesh = ops.tmin:0.05:ops.tmax;
plot(tmesh,Vtsig*100,'linewidth',2)
vline(0,'k')
xlabel(strcat(['time to ',alignlist_cpd{j}]))
ylabel(' CPD (%)')
title(strcat('CPD: ',groupnames_cpd{j}))
set(gca,'fontsize',15)



