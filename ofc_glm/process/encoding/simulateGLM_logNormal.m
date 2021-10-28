function [y,lambda,idx_struct,lambda_var] = simulateGLM_logNormal(c,c_cov,stim,ops,trial_id)
%generates synthetic data from a glm, but from lognormal dist
%output mean as from a log-normal distribution
%
%This is basically the core code from the gradient that calculates lambda
%but the spike times are kept in real time to capture the correct spike
%history filter response


%how many time-dependent covariates are causal/acausal
if ~isempty(stim.TDcausalvec)
    TDcaus_inds = find(stim.TDcausalvec(1,:)); %indices, causal filtering
    TDacaus_inds = find(stim.TDcausalvec(2,:)); %indices, acausal filtering
else
    TDcaus_inds = []; 
    TDacaus_inds = [];
end
ncaus = numel(TDcaus_inds); %# causal, time-dep. stimuli
nacaus = numel(TDacaus_inds); %# anticausal, time-dep. stimuli
nx = numel(stim.xmat); %# time-independent parameters, including bias

%stimuli. basisfun order: [spikes, stims]
lh = ops.lh; %spike history kernel length
bfh = stim.basefuns{1}; %spike filter basis functions
bfs = stim.basefuns{2}; %stimili filter basis functions. 
[dimh,~] = size(bfh); %dim, length of spike filter
%dims = zeros(ops.ncov-1,1);
dims = zeros(ops.nkern,1);
ls = zeros(ops.nkern,1);
%for j = 1:ops.ncov-1
for j = 1:ops.nkern
    [dims(j),ls(j)] = size(bfs{j}); %dim, length of stim filter
end
dim_stim = sum(dims(TDcaus_inds))+sum(dims(TDacaus_inds)); %total # basis functions for stimuli
%nc = dimh + dim_stim + nx; %total number of params

%% calculate the rate, lambda(t).TODO THIS SHOULD BE SIMULATE_GLM


%decide how spike timing and trial structure should be handled
trial_id_old = ops.Dat.trial_idx(trial_id); %full trial id, regardless of pseudotrials

if ops.usePseudotrial
    %spike and stimulus data is in continous time. find blocks of
    %continous time data
    [trial_idx] = parse_pseudotrials(ops,trial_id);    
    
else %trial-based structure for spike timing and stimuli
    %location of start and stop once vectorized. ensure they are col. vec.
    trial_idx = ops.Dat.trial_idx(trial_id);
    trs = reshape(ops.Dat.start_idx(trial_idx),numel(ops.Dat.start_idx(trial_idx)),[]);
    tre = reshape(ops.Dat.end_idx(trial_idx),numel(ops.Dat.end_idx(trial_idx)),[]);
    trial_idx = [trs,tre];
end

ntns = sum(diff(trial_idx,1,2)+1); %number of total time points
ns = size(trial_idx,1); %number of trials, or pseudotrials



%convolutions of stimuli. store conv with each basis function
stimh = zeros(dim_stim,ntns); 
%constant terms
constant_terms = zeros(ntns,1);
constant_terms_cov = zeros(ntns,1);
constant_params = c(end-nx+1:end);
%a fleshed out version of time-independent covariates. for gradient calc
xkron = zeros(nx,ntns);


%the causal convolutions for time-dependent stimuli
ind = 0;
%loop over all "trials"
for j = 1:ns
    nt = numel(trial_idx(j,1):trial_idx(j,2));
    
    
    ind_m = 0; %keeps track of where to place in stimh. iterates over basis fun
    for m = 1:ncaus
        bfs_m = bfs{m};
        stim_m = stim.TDstim{TDcaus_inds(m)};
        stim_mj = stim_m(trial_idx(j,1):trial_idx(j,2));
        
        for k = 1:dims(m)
            ind_m = ind_m+1;   
            stimh(ind_m,ind+1:ind+nt) = filter(bfs_m(k,:),1,stim_mj);
        end
    end
    
    %the anti-causal convolutions
    %ind_m = 0;
    for m = 1:nacaus
        bfs_m = bfs{m+ncaus};
        stim_m = stim.TDstim{TDacaus_inds(m)};
        stim_mj = stim_m(trial_idx(j,1):trial_idx(j,2));
        %ensure correct shape
        stim_mj = reshape(stim_mj,numel(stim_mj),[]);
       
        %stimulus, padded to stimulate event occuring ls steps in future
        ls_m = ls(m+ncaus);
        stim_padded = [stim_mj;zeros(ls_m,1)];   
        for k = 1:dims(m+ncaus)
            ind_m = ind_m + 1;
            sj = fliplr(bfs_m(k,:)); %time-reversed basis funcion
            sfilt = filter(sj,1,stim_padded);
            %stimh(ind_m,ind+1:ind+nt) = filter(sj,1,stim_padded);
            stimh(ind_m,ind+1:ind+nt) = sfilt(ls_m+1:end);
        end
    end
    
    %constant terms. this should have been treated as final rows
    %in the design matrix
    for m = 1:nx
        stim_m = stim.xmat{m};
        stim_m = stim_m(trial_idx(j,1):trial_idx(j,2));
        xkron(m,ind+1:ind+nt) = stim_m;
        constant_terms(ind+1:ind+nt) = constant_params(m)*stim_m;  
        constant_terms_cov(ind+1:ind+nt) = ...
            sum( (c_cov(end-nx+m,end-nx+m)*stim_m).*stim_m);
    end
       
    ind = ind + nt;
end

%constant params
stimparams = c(dimh+1:end-nx); %mean of kernel params
stimparams_cov = c_cov(dimh+1:end-nx,dimh+1:end-nx); %covariance of kernel params
stimconv = (stimparams*stimh)';
stimconv_cov = sum((stimparams_cov*stimh).*stimh)'; %new bit


%the actual simulation loop
%disp('simulating')
%nt = numel(stimconv);

%NOT beta tested
if ops.useSpike %use spike history kernel. proceed trial by trial
   
histparams = c(1:dimh);
histparams_cov =  c_cov(1:dimh,1:dimh);
histfilt = histparams*bfh; %history kernel
histfilt_cov = sum((histparams_cov*bfh).*bfh)';
    
lambda = zeros(ntns,1); %firing rate
lambda_var = zeros(ntns,1); %estimated variance of firing rate
y = zeros(ntns,1); %spike counts

ind = 0;
for m = 1:ns
    inds = (trial_idx(j,1):trial_idx(j,2));
    nt = numel(inds); %number time points of this trial
    trial_idx(m,:) = [ind+1,ind+nt];
    
    exparg = [constant_terms(inds)+stimconv(inds);zeros(lh,1)]; %mean of summed kernels
    exparg_cov = [constant_terms_cov(inds)+ stimconv_cov(inds);zeros(lh,1)]; %variance of summed kernels
    
    for j = 1:nt
        %lambda(ind+j) = exp(exparg(j));
        lambda(ind + j) = exp(exparg + exparg_cov/2);
        lambda_var(ind + j) = (exp(exparg_cov)-1).*exp(2*exparg+exparg_cov);
        
        
        y(ind+j) = poissrnd(lambda(ind+j)*ops.dt);
        if y(ind+j) > 0
            %add history kernel forward in time to exparg (the mean
            %component) and exparg_cov (the variance component)
            exparg(j:j+lh-1) = exparg(j:j+lh-1)+histfilt*y(ind+j);
            exparg_cov(j:j+lh-1) = exparg_cov(j:j+lh-1)+ histfilt_cov*y(ind+j)/2 ;
        end
    end
    ind = ind + nt;
end
   
    
else %no spike history kernel. don't simulate in time
    mu = stimconv + constant_terms;
    sigma_squared = stimconv_cov + constant_terms_cov;
    
    lambda = exp(mu + sigma_squared/2);
    %lambda = exp(stimconv+ constant_terms );
    y = poissrnd(lambda*ops.dt);   
    
    lambda_var = (exp(sigma_squared)-1).*exp(2*mu+sigma_squared);
    
end

%calulate start/end indices of original trials from trial_id in new output
%data frame.
%the important idx vectors in the new reference frame of simulated data
idx_struct = struct();
idx_struct.trial_idx = trial_id_old;
ns_old = numel(trial_id_old);

%start and end
idx_struct.start_idx = zeros(1,ns_old);
idx_struct.end_idx = zeros(1,ns_old);
ind = 0;
for j = 1:ns_old
    nt = numel(ops.Dat.start_idx(trial_id_old(j)):ops.Dat.end_idx(trial_id_old(j)));
    idx_struct.start_idx(j) = ind +1;
    idx_struct.end_idx(j) = ind + nt;
    ind = ind + nt;
end

%entercpoke, leavecpoke, choice
idx_struct.entercpoke_idx = idx_struct.start_idx - ...
    (ops.Dat.start_idx(trial_id_old) - ops.Dat.entercpoke_idx(trial_id_old));
idx_struct.leavecpoke_idx = idx_struct.start_idx - ...
    (ops.Dat.start_idx(trial_id_old) - ops.Dat.leavecpoke_idx(trial_id_old));
idx_struct.choice_idx = idx_struct.start_idx - ...
    (ops.Dat.start_idx(trial_id_old) - ops.Dat.choice_idx(trial_id_old));







