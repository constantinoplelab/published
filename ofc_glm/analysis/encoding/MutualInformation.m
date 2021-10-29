function [MI_td,shuffmean,shuffstd,shuffsig,Hs]...
        = MutualInformation(fname,ops)
%calculates mutual information between a stimulus or behavior of
%interest, and the firing rate, based on the underlying GLM
%
%%Inputs:
%   fname: glm fit file
%   ops: main structure with following info
%       masktype: 'test' or 'all' for held out of all data
%       condtype: string that is fet into condmask. for fig. title
%       condmask: function handle to return a conditional set of trials
%       alignment: 'start','leavecpoke','choice'
%       tmin: min time relative to start a trial
%       tmax: max time after start for a trial
%       nsamp: number of samples for building spike ct. distribution.
%
%
%Outputs,  nt time points in an event-aligned trial:
%   MI_td (1 x nt) matrix of model-predicted rate
%   shuffmean (nt x 1) vector of shuffled distribution mean
%   shuffstd (nt x 1) vector of shuffle distrbution standard deviation
%   Hs: Entropy of stimulus, in bits

f = load(fname); %load glm fit

tmesh = ops.tmin:f.ops.dt:ops.tmax;
nt = numel(tmesh);

nshuff = 500; %number shuffles for null distribution

%reference dist parts. 95% confidence interval boundary
shuffsig = zeros(nt,1);

%calculate firing rates according to log-normal distribution
%obtain timd-dependent mean and standard deviation of the log-rate, which
%follow a normal distribution
     
%decide on trial mask type: all of the data, or held-out test data
switch ops.masktype
    case 'test'
        trialmask = f.ops.test_ids;
    case 'all'
        trialmask = true(size(f.ops.test_ids));
    otherwise
        disp('no correct trial mask given ("test" or "all" . default to all data')
         trialmask = true(size(f.ops.test_ids));
end

%parse up trials and simulate from model
%log-normal one: lambda = exp(mu_t + sigma^2/2)
[~,lambda_ln,~] = cleantrials(f,trialmask,ops.alignment, ops.tmin,ops.tmax,struct('simtype','lognormal'));
%point-estimate normal one: lambda = exp(mu_t)
[~,lambda_glm,trial_idx_clip] = cleantrials(f,trialmask,ops.alignment, ops.tmin,ops.tmax,[]);


%filter out some trials so only relevant conditions are present
%this is a way to get at mutual information at a finer level
if isfield(ops,'condmask')             
    %find set of all relevant trials
    condmask_sum = zeros(size(ops.condmask{1}));
    for k = 1:numel(ops.condmask)
        condmask_sum = condmask_sum + ops.condmask{k};
    end
    condmask_sum = boolean(condmask_sum);     
    condtrial = f.ops.Dat.trial_idx(condmask_sum);
    condmask_all = ismember(trial_idx_clip,condtrial);  

    %remove irrelevant trials from the original condmasks
    for k = 1:numel(ops.condmask)
        ops.condmask{k} = ops.condmask{k}(condmask_all);
    end

else %no conditional masking
    condmask_all = true(size(trial_idx_clip));
end

lambda_ln = lambda_ln(condmask_all,:);
lambda_glm = lambda_glm(condmask_all,:);
%dat = dat(condmask_all,:);

%basics
maxct = 10; %max spike number per bin allowed (200 Hz)
n_cond = numel(ops.condmask);
[ns,nt] = size(lambda_ln);
pvec = zeros(n_cond,1);

%% the sampling approach

shuff_trials = zeros(nshuff,ns);
for j = 1:nshuff
    shuff_trials(j,:) = randperm(ns);
end


%extract mu,sigma^2;
mu_t = log(lambda_glm);
sigma_sq_t = 2*log(lambda_ln./(lambda_glm)); 

%sample lambda values by particle filter or normal throguh exponential
nsamp = ops.nsamp;
lambda_samp = zeros(nsamp,size(mu_t,1),size(mu_t,2));

lamj =  exp(mvnrnd(mu_t(:),sigma_sq_t(:)',nsamp));

for j = 1:nsamp
    lambda_samp(j,:,:) = reshape(lamj(j,:),ns,nt);
end


%propagate to poisson, and build count distribution
%counts with rate parameter drawn from above samples of lambda
y_samps = zeros(nsamp,ns,nt);
for j = 1:ns
    for k = 1:nt
        y_samps(:,j,k) = poissrnd(lambda_samp(:,j,k)*f.ops.dt);
    end
end


%create the conditional distributions over spike counts
y_counts = zeros(maxct,ns,nt);
for j = 1:ns
    for k = 1:nt
        y_counts(:,j,k) = histcounts(y_samps(:,j,k),0:maxct);
    end
end

Hj = zeros(n_cond,nt);
Hj_shuff = zeros(nshuff,n_cond,nt);



for j = 1:n_cond
    %mask out each condition as categories of the stimulus distribution
    %account for edge clipping
    condtrial = f.ops.Dat.trial_idx(ops.condmask{j});
    condmask = ismember(trial_idx_clip,condtrial);
    %stimulus probability
    pvec(j) = sum(condmask)/ns;

    %build conditional distribution
    p_j = squeeze(sum(y_counts(:,condmask,:),2));
    ynorm = sum(p_j,1);
    p_j = p_j./ynorm;

    %take actual expectation of -log(p), avoiding zeros that make nans
    Hj(j,:) = nansum(-p_j.*log(p_j));

    %do the same for shuffled distribution       
    for k = 1:nshuff
        trial_idx_clip_shuff = trial_idx_clip(shuff_trials(k,:));           
        condtrial = f.ops.Dat.trial_idx(ops.condmask{j});
        condmask = ismember(trial_idx_clip_shuff,condtrial);

        %build conditional distribution
        p_jshuff = squeeze(sum(y_counts(:,condmask,:),2));
        ynorm = sum(p_jshuff,1);
        p_jshuff = p_jshuff./ynorm;

        %take actual expectation of -log(p), avoiding zeros that make nans
        Hj_shuff(k,j,:) = nansum(-p_jshuff.*log(p_jshuff));                      
    end       

end

MI_shuff_samp = zeros(nshuff,nt);

%what will H be? via the proper marginal
p_marg = squeeze(sum(y_counts,2));
p_marg_norm = sum(p_marg,1);
p_marg = p_marg./p_marg_norm;
H_samp = -nansum(p_marg.*log(p_marg));
Hcond_samp = pvec'*Hj;
MI_td = H_samp- Hcond_samp;
MI_td = MI_td/log(2);

%entropy of stimulus
Hs = -pvec'*log(pvec)/log(2);

%for shuffled dist
for k = 1:nshuff
    Hcond_samp_shuff = pvec'*squeeze(Hj_shuff(k,:,:));
    MI_td_shuff = H_samp- Hcond_samp_shuff;
    MI_shuff_samp(k,:) = MI_td_shuff/log(2);  
end

shuffmean = nanmean(MI_shuff_samp);
shuffstd = nanstd(MI_shuff_samp);

%for every time point, calculate the 95% CI
for l = 1:nt
    mh = MI_shuff_samp(~isnan(MI_shuff_samp(:,l)),l);
    if ~isempty(mh)
        [h,x] = histcounts(mh,20,'Normalization','cdf');
        g95 = find(h > 0.975);
        shuffsig(l) = x(g95(1)); %first instance
    end
end


    










