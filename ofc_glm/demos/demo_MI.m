%mutual information between spikes and trial-level behavior/outcome stimuli

%% choose a neuron and its glm fit, as well as place to save


%% paths to stuff on local machine
datadir = '~/projects/ofc/data/published/'; %where raw data is stored


%where will feature space data be stored. and where glm fit file is stored
%can change them to be separate directories if you would like
savedir = '/Users/dhocker/projects/ofc/results/demo/'; 

codedir = '~/projects/constantinoplelab/published/ofc_glm/'; %where code lives

file_idx = 170; %neuron to compute MI on

addpath(genpath(codedir));

%this is the name of the output glm fit from the demo, but can in general
%gen anything
fname_glm = strcat(savedir,'fit_stim2020903_N',num2str(file_idx),'.mat');


%% load the A struct of concatenated data 
%to help with creating trial contingency masks

D = load(strcat(datadir,'concatdata_ofc_noOptOut.mat'));
A = D.A;
clear D

%% helper cells to define which kernels and trial contingencies to use 

%groupings for conditional PSTHs and mutual information
condgroups = {{'prevWin','prevLoss','prevOptOut'},...
            {'left','right'},...
            {'win','loss'},...
            {'risky','safe'}};     
aligngroups =   {'start','leavecpoke','choice','choice'};
groupnames = {'rewardHistory','LR','reward','riskySafe'};

%generate trial masks for each type of trial contengincy
condmask = cell(numel(condgroups),1);
for j = 1:numel(condgroups)    
    %the common case, sets of binary conditions, one per trial, and easily
    %distinguished as a case from condmaskLibrary with a string handle
    if iscell(condgroups{j})
        cj = cell(numel(condgroups{j}),1);
        for k = 1:numel(condgroups{j})
            cj{k} = condmaskLibrary(condgroups{j}{k} ,file_idx,A,struct('removeOptOut',true));
        end
    else %case for scalar condition with single 'string' handle. e.g., 'vol'
        cj = condmaskLibrary(condgroups{j} ,file_idx,A,struct('removeOptOut',true));
    end
    
    condmask{j} = cj;
end

%time ranges for each alignment type, in seconds from alignment event
tminvec = [-2,-4,-4];
tmaxvec = [4,4,4];
alcell = {'start','leavecpoke','choice'};

%% Mutual information.


%this requires groups of condtypes and condmasks that deal with the same
%situation


j = 1; %index in the helper lists above. which cpd calc to perform


%should models be a cell-struct to plot multiple? yes. for now just one


ops = struct();
ops.masktype = 'all'; %test data or all data?
ops.condtype = condgroups{j};
ops.condmask = condmask{j}; %contingency masks for trial types
ops.alignment = aligngroups{j}; %what MI will be aligned to   
ops.tmin = tminvec(strcmp(aligngroups{j},alcell));%trial start time (s)
ops.tmax = tmaxvec(strcmp(aligngroups{j},alcell));  %trial end time (s)
ops.simtype = 'glm'; %how to simulate form GLM
ops.nsamp = 500; %number of samples from Poiss. dist in each time bin


disp('beginning MI calc')
[MI,shuffmean,shuffstd,shuffsig,Hs] = MutualInformation(fname_glm,ops);

%% this was the raw mutual information being calculated. significance correct
MI_sig = MI;
MI_sig(MI_sig < shuffsig) = 0;

    
%% save   
savename = strcat([savedir,'MutualInformation_N',num2str(file_idx),'_',groupnames{j},'.mat']);
save(savename,'MI_sig','MI','shuffmean','shuffstd','Hs',...
    'shuffsig','ops')

%% plot
dt=0.05;
tmesh = ops.tmin:dt:ops.tmax;
figure(89)
clf
plot(tmesh,MI_sig/dt,'linewidth',2);
vline(0,'k')
xlabel(strcat(['time to ',ops.alignment]))
ylabel(' MI (bits/s)')
title(strcat('MI: ',groupnames{j}))
set(gca,'fontsize',15)

    