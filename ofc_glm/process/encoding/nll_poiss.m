function [f,G,H,evidence,penalty,f_offset] = nll_poiss(st,c,stim,ops,trial_id)
%calcs neg. log-likelihood and gradient, hessian
%for glm with parametrized spike history filter and multiple filtered stimuli
%
%
%Inputs:
%st: either i) [nt,ns] matrix of spike times, nt=# time points,ns=# samples
%    or     ii) [nt x ns,1] vector of spike times. a continuous time form
%c: [nc,1] control parameters
%stim: struct for covariates to model. see parse_stimuli_* for more info
%ops: main structure for control flow arguments. see glm_defaultOps.m for
%all options
%trial_id: [ns,1] boolean of which trials to use for calculation
%
%Outputs
%f: negative log-likelihood(or posterior)
%G: [nc,1] graident of negative log-likelihood(or poseterior) vector w.r.t 
    %control params in c
%H: [nc,nc] Hessian matrix
%evidence: log of evidence. useful for hyperparameter optimization

%how many time-dependent covariates are causal/acausal
if ~isempty(stim.TDcausalvec)
    TDcaus_inds = find(stim.TDcausalvec(1,:)); %indices, causal filtering
    TDacaus_inds = find(stim.TDcausalvec(2,:)); %indices, acausal filtering
else
    TDcaus_inds = []; 
    TDacaus_inds = [];
end
nTD_caus = numel(TDcaus_inds); %# causal, time-dep. stimuli
nTD_acaus = numel(TDacaus_inds); %# anticausal, time-dep. stimuli
nx = numel(stim.xmat); %# time-independent parameters, including bias

assert(nTD_caus+nTD_acaus==ops.nkern,'ops.nkern set incorrectly');

%stimuli. basisfun order: [spikes, stims]
bfh = stim.basefuns{1}; %spike filter basis functions
bfs = stim.basefuns{2}; %stimili filter basis functions. 
[dimh,~] = size(bfh); %dim, length of spike filter

%TODO: CHANGE TO NUMBER OF KERNELS
%dims = zeros(ops.ncov-1,1);
%ls = zeros(ops.ncov-1,1);
dims = zeros(ops.nkern,1);
ls = zeros(ops.nkern,1);
for j = 1:ops.nkern
    [dims(j),ls(j)] = size(bfs{j}); %dim, length of stim filter
end
dim_stim = sum(dims(TDcaus_inds))+sum(dims(TDacaus_inds)); %total # basis functions for stimuli
nc = dimh + dim_stim + nx; %total number of params

%% calculate the rate, lambda(t).TODO THIS SHOULD BE SIMULATE_GLM


%decide how spike timing and trial structure should be handled
if ops.usePseudotrial
    %spike and stimulus data is in continous time. find blocks of
    %continous time data
    [trial_idx] = parse_pseudotrials(ops,trial_id);
    st_all = st;
   
else %trial-based structure for spike timing and stimuli
    %location of start and stop once vectorized. ensure they are col. vec.
    trial_idx = ops.Dat.trial_idx(trial_id);
    trs = reshape(ops.Dat.start_idx(trial_idx),numel(ops.Dat.start_idx(trial_idx)),[]);
    tre = reshape(ops.Dat.end_idx(trial_idx),numel(ops.Dat.end_idx(trial_idx)),[]);
    trial_idx = [trs,tre];
    
    %keep copy of original spike time for correct indexing
    st_all = reshape(st,numel(st),[]); %reshape into vectorized form
  
end

ntns = sum(diff(trial_idx,1,2)+1); %number of total time points
ns = size(trial_idx,1); %number of trials, or pseudotrials  
st = zeros(ntns,1); %only trial-relevant data for LL, grad, etc.
ind = 0;
for j = 1:ns
   nt = numel(trial_idx(j,1):trial_idx(j,2));
   st(ind + 1: ind + nt) = st_all(trial_idx(j,1):trial_idx(j,2));
   ind = ind + nt;
end

%using autoregressive spike history kernel? convolve spike train
if ops.useSpike
    %convolute spike filter with spikes
    %perform filtering rowwise, then reshape to avoid trials interacting with
    histparams = c(1:dimh);
    histfilt = histparams*bfh;
    
    %convoled spike train with spike history kernel
    histconv = zeros(ntns,1);
    %convolved spike train with each basis function element
    yh = zeros(dimh,ntns);
    
    ind = 0;
    for j = 1:ns
        st_j = st_all(trial_idx(j,1):trial_idx(j,2));
        nt = numel(st_j); %# time points for this data
        histconv(ind+1:ind+nt) = filter(histfilt,1,st_j);

        for k = 1:dimh
        yh(k,ind+1:ind+nt) = filter(bfh(k,:),1,st_j);
        end
        
        ind = ind + nt;
    end
         
else
    yh = zeros(dimh,ntns);
    histconv = zeros(ntns,1);
end       


%convolutions of stimuli. store conv with each basis function
stimh = zeros(dim_stim,ntns); 
%constant terms
constant_terms = zeros(ntns,1);
constant_params = c(end-nx+1:end);
%a fleshed out version of time-independent covariates. for gradient calc
xkron = zeros(nx,ntns);


%the causal convolutions for time-dependent stimuli
ind = 0;
%loop over all "trials"
for j = 1:ns
    nt = numel(trial_idx(j,1):trial_idx(j,2));
    
    ind_m = 0; %keeps track of where to place in stimh. iterates over basis fun
    for m = 1:nTD_caus
        %TODO: CHANGE to simpler indexing now that all kernels have a bf
        %bfs_m = bfs{TDcaus_inds(m)};
        bfs_m = bfs{m};
        stim_m = stim.TDstim{TDcaus_inds(m)};
        stim_mj = stim_m(trial_idx(j,1):trial_idx(j,2));
        
        %for k = 1:dims(TDcaus_inds(m))
        for k = 1:dims(m)
            ind_m = ind_m+1;   
            stimh(ind_m,ind+1:ind+nt) = filter(bfs_m(k,:),1,stim_mj);
        end
    end

    %the anti-causal convolutions
    %ind_m = 0;
    for m = 1:nTD_acaus
        %TODO: CHANGE to simpler indexing now that all kernels have bf
        %bfs_m = bfs{TDacaus_inds(m)};
        bfs_m = bfs{m+nTD_caus};
        stim_m = stim.TDstim{TDacaus_inds(m)};
        stim_mj = stim_m(trial_idx(j,1):trial_idx(j,2));
        %ensure correct shape
        stim_mj = reshape(stim_mj,numel(stim_mj),[]);
       
        %stimulus, padded to stimulate event occuring ls steps in future
        %ls_m = ls(TDacaus_inds(m));
        ls_m = ls(m + nTD_caus);
        stim_padded = [stim_mj;zeros(ls_m,1)];   
        %for k = 1:dims(TDacaus_inds(m))
        for k = 1:dims(nTD_caus+m)
            ind_m = ind_m + 1;
            sj = fliplr(bfs_m(k,:)); %time-reversed basis funcion
            sfilt = filter(sj,1,stim_padded);
            %stimh(ind_m,ind+1:ind+nt) = filter(sj,1,stim_padded);
            stimh(ind_m,ind+1:ind+nt) = sfilt(ls_m+1:end);
        end
    end
    
    %constant terms
    for m = 1:nx
        stim_m = stim.xmat{m};
        stim_m = stim_m(trial_idx(j,1):trial_idx(j,2));
        xkron(m,ind+1:ind+nt) = stim_m;
        constant_terms(ind+1:ind+nt) = constant_params(m)*stim_m;        
    end
       
    ind = ind + nt;
end

%constant params
stimparams = c(dimh+1:end-nx);
stimconv = (stimparams*stimh)';


lambda = exp(histconv+stimconv+ constant_terms )*ops.dt;


%% calculate terms for use in gradient and hess

%regularizer penalty for L2
if isfield(ops,'gamma'), gamma = ops.gamma; else, gamma = 0; end
%regularizer penalty for L1
if isfield(ops,'beta'), beta = ops.beta; else, beta = 0; end
%regularizer = ops.regularizer; %not importnat yet

nzind = find(lambda); %nonzero indices. avoid -inf


f = nansum(lambda(nzind))-st(nzind)'*log(lambda(nzind));
f_offset = nansum(lambda(nzind))-st(nzind)'*log(lambda(nzind)) - sum(st(nzind))*log(1/ops.dt);

switch ops.regularizer
    case 'L1'
        penalty = sum(beta*abs(c));
    case 'L2'
        penalty = sum(gamma/2*c.^2);
    case 'L1+L2'
        penalty = [sum(gamma/2*c.^2), sum(beta*abs(c))];
end

%penalty = sum(gamma/2*c(1:end-1).^2) + sum(beta*abs(c(1:end-1))); %allow bg rate to be un-penalized
%f = f/(ns*nt);

f = f+ sum(penalty);


%order: spike history, stim filter, constant functions 
%xkron =  kron(ones(1,nt),stim.xmat); %size nx, ns*nt
G = [yh*(lambda-st); stimh*(lambda-st); ...
    xkron*(lambda-st)];

G = G + gamma*c' +  beta*(c'~= 0);

%G = G + gamma*[c(1:end-1)';0] +  beta*[(c(1:end-1)'~= 0);0]; not reg. on
%last term
%G = G/(ns*nt);

%hessian
%param suffixes: y:spike-filter, stim:stim-filter, x:constant-param

%smaller memoery requirements
%block diagonal terms
kron_lambda_sp = kron(lambda,ones(1,dimh));
kron_lambda_nx = kron(lambda,ones(1,nx));

Hy = yh*(kron_lambda_sp.*yh');
Hstim = stimh*(kron(lambda,ones(1,dim_stim)).*stimh'); 
Hx = xkron*(kron_lambda_nx.*xkron');

%block off-digonals        
Hyx = yh*(kron_lambda_nx.*xkron'); %size dimh,nx 
Hstimx = stimh*(kron_lambda_nx.*xkron'); %size dim_stim,nx 
Hstimy = stimh*(kron_lambda_sp.*yh'); %size dim_stim, dimh

H = zeros(nc);      
%starting indices
i1 = 0;
i2 = dimh;
i3 = dimh+dim_stim;
%digonal terms
H(i1+1:i1+dimh ,   i1+1:i1+dimh)        = Hy;
H(i2+1:i2+dim_stim , i2+1:i2+dim_stim)  = Hstim;

H(i3+1:i3+nx ,       i3+1:i3+nx)        = Hx;
%H(i3+1:i3+nx ,       i3+1:i3+nx)        = 0; %allow bg rate as unpenalized

%off-digonal
H(i1+1:i1+dimh ,   i3+1:i3+nx)          = Hyx;
H(i3+1:i3+nx ,       i1+1:i1+dimh)      = Hyx';
H(i2+1:i2+dim_stim , i1+1:i1+dimh)      = Hstimy;
H(i1+1:i1+dimh ,   i2+1:i2+dim_stim)    = Hstimy';
H(i2+1:i2+dim_stim ,  i3+1:i3+nx)       = Hstimx;
H(i3+1:i3+nx ,       i2+1:i2+dim_stim)  = Hstimx';

%H = H/(ns*nt);
switch(ops.regularizer)
    case 'L2'
        H = H + gamma*eye(size(H)); 
    case 'L1+L2'
        H = H + gamma*eye(size(H));         
end

%old way
%evidence = -f - 0.5*log(det(H))+ 0.5*nc*log(2*pi); 
%new way
switch(ops.regularizer)
    case 'L1'
        %evidence = -f+nc*log(beta/2)-0.5*log(det(H));

        evidence = -f-0.5*log(det(H));
    case 'L2'
        %i think 0.5*log(gamma) is a typo, since this is already
        %is contained in H from above.
        %evidence = -f+0.5*log(gamma)-0.5*log(det(H));

        %general correct form?
        evidence = -f-0.5*log(det(H));
    otherwise
        evidence = -f-0.5*log(det(H));
end      

        


