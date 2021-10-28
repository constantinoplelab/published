function [Vt,Vtshuff,shuffmean,shuffstd] = CPD_noOpt(file_idx,ops)
%a non-optimization-based form of coefficient of partial determination
%calculation, a way to assess relative changes in %variance explained due
%to a particular kernel in a model
%Inputs:
%   file_idx: neuron index in concatenated A data structure, or index in
%           fnames
%   ops: main structure with following info
%       fitdir: where fit files live
%       datadir: where raw data files live
%       namebase: 'prefix of specific fit file. often is
%           "fit_stim{STIMTYPE}N"
%       masktype: 'test' or 'all' for held out of all data
%       alignment: cell of 'start','leavecpoke','choice'. one for each
%               kernel
%       shufftype: 'behavior' to combine several covariates, or 'omit' to
%                   remove them outright
%       cpdname: the name of the CPD being tests
%       cpdlist: list of kernel names to be either omitted or combined
%
%   makeFig: boolean to turn on (true) or off (false) plotting
%
%Outputs, for ns trial, nt time points
%   Vt: (nt x 1) CPD
%   stimlegend: (nshuff x nt) matrix of shuffled CPD values


fname = strcat([ops.fitdir,ops.namebase,num2str(file_idx),'.mat']); 
f = load(fname);
stim = f.stim;
alignment = ops.alignment;

%change path of load directory to make sure it points to correct data
tmpstr = split(f.ops.dataname,'/');
f.ops.dataname = strcat([ops.datadir,tmpstr{end}]);  


%which data to use in calculation
switch ops.masktype
    case 'test'
        trialmask = f.ops.test_ids;
    case 'all'
        trialmask = true(size(f.ops.test_ids));
end

TDcaus_inds = find(f.stim.TDcausalvec(1,:));
TDacaus_inds = find(f.stim.TDcausalvec(2,:));
TD_inds = [TDcaus_inds,TDacaus_inds];

%which kernels will be combined/omitted
stiminds = [];
for j = 1:numel(ops.cpdgroup)
    stiminds = cat(1,stiminds,TD_inds(strcmp(ops.cpdgroup{j},f.stim.stimlegend)));
end
stiminds = unique(stiminds);
ngroup = numel(stiminds); %this handles correct case when both causal and acausal are affected by a stimulus


if isempty(stiminds) %do not do calcualation. not present
    disp(strcat(['kernels for ',ops.cpdname,' not present in model ',ops.modelnames]))
    Vt = [];
    Vtshuff = [];
    shuffmean = [];
    shuffstd = [];
    return;
end

%stimulate
disp('simulate full')
[rate_dat,lambda_allkern,trial_idx_clip] = cleantrials(f,trialmask, alignment,ops.tmin,ops.tmax,ops);

%% which trials belong to which condition?

if isfield(ops,'condmask')               
    if iscell(ops.condmask)
        condmask = zeros(numel(trial_idx_clip),numel(ops.condmask));
        for k = 1:numel(ops.condmask) 
            condtrial = f.ops.Dat.trial_idx(ops.condmask{k});  
            condmask(:,k) = ismember(trial_idx_clip,condtrial);
        end
    else
        condmask = zeros(numel(trial_idx_clip),1);
        condtrial = f.ops.Dat.trial_idx(ops.condmask);  
        condmask(:) = ismember(trial_idx_clip,condtrial);
    end
else %no conditional masking
    %condmask_all = true(size(trial_idx_clip));
    condmask = true(size(trial_idx_clip));
end


%% decide how to omit the kernels of interest to check behavior
disp('leaveout kernels')
fwo = f;
switch ops.shufftype
    case 'behavior'
        fwo.ops.omitBehavior = ops.cpdgroup;
        for j = 1:size(fwo.ops.omitBehavior,1)
            fwo.ops.omitBehavior{j,end+1} = strcat(['nb ',num2str(j)]);
        end
    case 'omit'
        fwo.ops.omitKernels = ops.cpdgroup;
end


%TEMP REMOVE: updates ops sizing if a special form of basis has been used
%if strcmp(f.ops.basistype,'ofc_adaptive_plus')
%    fwo.ops.dimst(:) = 9;   
%end

[Dat_wo,stim_wo] = parse_stimuli_wrapper(fwo.ops);
[basefun,basefunst] = chooseBasis_ofc(fwo.ops);
stim_wo.basefuns = {basefun,basefunst};

%TEMP REMOVE: updates ops sizing if a special form of basis has been used
%if strcmp(f.ops.basistype,'ofc_adaptive_plus')
%    nbnew = size(stim_wo.basefuns{2}{1},1);
%    fwo.ops.dimst(:) = nbnew;   
%end





%find indices of kernels being omitted/merged
stiminds_vec = zeros(numel(ops.cpdgroup),1); %where in wML are you?
avgroups_vec = zeros(numel(ops.cpdgroup),1); %are you causal(1) or acausal(2)? for averaging
for j = 1:numel(ops.cpdgroup)
    [stimind,~,covops] = stimulusInfo(stim,(ops.cpdgroup{j}));
    stiminds_vec(j) = stimind;
    if covops.isCausal
        avgroups_vec(j) = 1;
    else
        avgroups_vec(j) = 2;
    end
end

switch ops.shufftype
    case 'behavior'
        %ASSUMES SAME BASIS FUNCTIONS AND
        %KERNELS. no other way to average right now
        
        %find and average params for causal, acausal
        stiminds_causal = stiminds_vec(avgroups_vec==1);
        paraminds_causal = [];
        for j = 1:numel(stiminds_causal)
            ind = sum(f.ops.dimst(1:stiminds_causal(j)-1));
            paraminds_causal = cat(1,paraminds_causal,ind+1:ind+f.ops.dimst(1:stiminds_causal(j)));
        end
        
        stiminds_acausal = stiminds_vec(avgroups_vec==2);
        paraminds_acausal = [];
        for j = 1:numel(stiminds_acausal)
            ind = sum(f.ops.dimst(1:stiminds_acausal(j)-1));
            paraminds_acausal = cat(1,paraminds_acausal,ind+1:ind+f.ops.dimst(1:stiminds_acausal(j)));
        end  
        
        
        %to the param updating for each fold
        for j = 1:numel(f.wML)
            
            %average causal ones
            wMLcaus = mean(fwo.wML{j}(paraminds_causal));
            wMLacaus = mean(fwo.wML{j}(paraminds_acausal));
            
            sizeops = struct(); %use size info of first stim
            sizeops.ls = f.ops.ls(stiminds_vec(1));
            sizeops.dimst = f.ops.dimst(stiminds_vec(1));
            sizeops.dimbcst = f.ops.dimbcst(stiminds_vec(1));
            
            %remove all old covariates
            wMLnew = fwo.wML{j};
            opsnew = fwo.ops;
            remove_inds = sort(stiminds_vec);
            for m = 1:numel(remove_inds) %loops over every covariate except last one added (the nonbeh one)
                [wMLnew,opsnew] = removeParams(wMLnew,remove_inds(m),opsnew);
                remove_inds = remove_inds-1; %index will have shifted by one since a covariate removed
            end
            
            %put the new ones in
            if ~isnan(sum(wMLcaus))
                [stimind,~,~] = stimulusInfo(stim_wo,fwo.ops.omitBehavior{end});
                [wMLnew,opsnew] = addParams(wMLnew,wMLcaus,stimind,opsnew,sizeops);
            end
            
            if ~isnan(sum(wMLacaus))
                %is this the first covarate or second?
                if ~isnan(sum(wMLcaus))
                    [stimind,~,~] = stimulusInfo(stim_wo,strcat(fwo.ops.omitBehavior{end},' (ac)'));
                else
                    [stimind,~,~] = stimulusInfo(stim_wo,fwo.ops.omitBehavior{end});
                end
                [wMLnew,opsnew] = addParams(wMLnew,wMLacaus,stimind,opsnew,sizeops);
            end
            
            fwo.wML{j} = wMLnew;
            
        end
        fwo.ops = opsnew;
        

        
    case 'omit' %collect and keep all kernels not related to stiminds
        
        
        %remove all old covariates
        for j = 1:numel(f.wML)
            wMLnew = fwo.wML{j};
            opsnew = fwo.ops;
            remove_inds = sort(stiminds_vec);
            for m = 1:numel(remove_inds) %loops over every covariate except last one added (the nonbeh one)
                [wMLnew,opsnew] = removeParams(wMLnew,remove_inds(m),opsnew);
                remove_inds = remove_inds-1; %index will have shifted by one since a covariate removed
            end
            
            fwo.wML{j} = wMLnew;
            
        end
        
        fwo.ops = opsnew;
                    
end

fwo.ops.Dat = Dat_wo;
fwo.stim = stim_wo;
fwo.ops.stim = stim_wo; %override and use this stim for cleantrials


disp('simulate reduced')
%simulate
[~,lambda_allkern_nobeh] = cleantrials(fwo,trialmask, alignment,ops.tmin,ops.tmax,ops);

%% calulate true cpd

%scale rates and predictions by relative number 
ntrial = sum(condmask,'all'); %how many total trials
ncond = size(condmask,2);
for k = 1:ncond
    ck = condmask(:,k);
    %wk = sum(ck)/ntrial; %weight based on proporition of prevelance
    wk = ntrial/(sum(ck)*ncond); %weight to balance across conditions
    ck = boolean(ck);
    rate_dat(ck,:) = wk*rate_dat(ck,:);
    lambda_allkern_nobeh(ck,:) = wk*lambda_allkern_nobeh(ck,:);
    lambda_allkern(ck,:) = wk*lambda_allkern(ck,:);
end


SSE_reduced_t = sum((rate_dat-lambda_allkern_nobeh).^2);
SSE_full_t = sum((rate_dat-lambda_allkern).^2);

Vt =  (SSE_reduced_t-SSE_full_t)./(SSE_reduced_t);

%% handle the shuffle. do first iteration outside loop with hyperparam opt
disp('begin shuffle')
%shuffled version 
%nshuff = 100;
nshuff = 500;

%shuffled trial indices.nzinds is all time indices where relevant stim were
%doing things
rng('default');
nz_inds = [];
for j = 1:numel(stiminds)
    nz_inds = cat(1,nz_inds,find(stim.TDstim{stiminds(j)}));
end
nz_inds = unique(nz_inds); 
nt_shuff2 = numel(stim.TDstim{1});
shmat = zeros(nshuff,ngroup,nt_shuff2); %scramble across two covariates
shmat2 = zeros(nshuff,nt_shuff2); %scrambe timing of individual covariate

for m = 1:nshuff
    for k = 1:numel(nz_inds)
        shmat(m,:,nz_inds(k)) = randperm(ngroup);
    end
       
    shmat2(m,:) = randperm(nt_shuff2,1); %needed?
end

tmesh = ops.tmin:f.ops.dt:ops.tmax;
nt = numel(tmesh);
Vtshuff = zeros(nshuff,nt);


%% shuffled version, do the shuffles with that optimized hyperparameter
disp('here')
for m = 1:nshuff
    disp(m)
    fshuffm = f;
    %shuffle
    stim_shuff = stim;
    TD = stim_shuff.TDstim(stiminds);
    TDnew = TD;
    switch ops.shufftype
        case 'behavior'                   
            for j = 1:ngroup
                for k = 1:ngroup
                    mask = shmat(m,j,:)==k;
                    TDnew{j}(mask) =TD{k}(mask);
                end
            end
        case 'omit'
            for j = 1:ngroup
                TDnew{j} = TDnew{j}(shmat2(m,:));
            end
    end
    stim_shuff.TDstim(stiminds) = TDnew;    
    
    fshuffm.stim = stim_shuff;
    fshuffm.ops.stim = stim_shuff; %overwrites things to use this stim


    %simulate
    [~,lambda_shuff] = cleantrials(fshuffm,trialmask, alignment,ops.tmin,ops.tmax,fshuffm.ops);
    
    for k = 1:ncond
        ck = condmask(:,k);
        %wk = sum(ck)/ntrial; %weight based on proporition of prevelance
        wk = ntrial/(sum(ck)*ncond); %weight to balance across conditionsg
        ck = boolean(ck);
        lambda_shuff(ck,:) = wk*lambda_shuff(ck,:);
    end
    
    %CPD
    SSE_shuff = sum((rate_dat-lambda_shuff).^2);
    Vtshuff(m,:) =  (SSE_reduced_t-SSE_shuff)./(SSE_reduced_t);

    
end

shuffmean = nanmean(Vtshuff);
shuffstd = nanstd(Vtshuff);








   
    




    