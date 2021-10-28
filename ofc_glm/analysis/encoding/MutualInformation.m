function [MI_tdList,shuffmean,shuffstd,MI_ub_tdList,...
        shuffmean_ub, shuffstd_ub,shuffsig,shuffsig_ub,Hs, tmesh,fig,...
        MI_glm_ub,MI_emp] = MutualInformation(file_idx,ops,makeFig)
%calculates mutual information between a stimulus or behavior of
%interest, and the firing rate, based on the underlying GLM
%
%%Inputs:
%   file_idx: neuron index in concatenated A data structure, or index in
%           fnames
%   ops: main structure with following info
%       fitdir: where fit files live
%       datadir: where raw data files live
%       namebase: 'prefix of specific fit file. often is
%           "fit_stim{STIMTYPE}N"
%       modelnames: readable form of namebase, for labels
%       masktype: 'test' or 'all' for held out of all data
%       condtype: string that is fet into condmask. for fig. title
%       condmask: function handle to return a conditional set of trials
%       alignment: 'start','leavecpoke','choice'
%       tmin: min time relative to start a trial
%       tmax: max time after start for a trial
%       fignum: for numbering figures
%       nsamp: number of samples for building spike ct. distribution. 1000?
%       stimname: name for covariate/behavior being analyzed
%
%   makeFig: boolean to turn on (true) or off (false) plotting
%
%Outputs,  nt time points in an event-aligned trial:
%   MI_tdlist (nmod x nt) matrix of model-predicted rate
%   shuffmean (nt x 1) vector of shuffled distribution mean
%   shuffstd (nt x 1) vector of shuffle distrbution standard deviation
%   Hs: Entropy of stimulus, in bits
%   tmesh: (nt x 1) timing of trial
%   fig: output figure handle
%also lower limits via an approximation with Jensen's inequality
%note: the shuffled mean is only saved for the last model in the list
%the appear to have very similar sizes when investigating a few simple
%cells


if iscell(ops.namebase)
    n_mods = numel(ops.namebase);
    flist = cell(n_mods,1);
    for k = 1:n_mods
        fname = strcat([ops.fitdir,ops.namebase{k},num2str(file_idx),'.mat']); 
        flist{k} = load(fname);
    end
    %colormat = distinguishable_colors(n_mods);
    colormat = linspecer(n_mods,'qualitative'); %color schemes
else %please don't put just one inside the cell
    n_mods = 1;
    fname = strcat([ops.fitdir,ops.namebase,num2str(file_idx),'.mat']); 
    flist = {load(fname)};
    colormat = [1,0,0];
end

%make figure
if makeFig
    if isfield(ops,'fignum')
        fig = figure(ops.fignum);
        clf
        hold on
    else
        fig = figure;
        clf
    end
end

%structure to save everything
tmesh = ops.tmin:flist{1}.ops.dt:ops.tmax;
nt = numel(tmesh);

nshuff = 500; %number shuffles for null distribution


MI_tdList = zeros(n_mods,nt);
MI_ub_tdList = zeros(n_mods,nt); %upper bound version


shuffmean = zeros(nt,1);
shuffstd = zeros(nt,1);
shuffmean_ub = zeros(nt,1); %upper bound version
shuffstd_ub = zeros(nt,1);

%95% confidence interval boundary
shuffsig = zeros(nt,1); %for sampling
shuffsig_ub = zeros(nt,1); %for upper bound
shuffsig_ub_glm = zeros(nt,1); %for glm upper bound

%for the upper bound based on a point estimate GLM model.
%might be formally more correct, since inference was done to 
%optimize the point estimate
MI_ub_td_glm_List = zeros(n_mods,nt);
shuffmean_ub_glm = zeros(nt,1); %upper bound version
shuffstd_ub_glm = zeros(nt,1);

%emperical distribution form
MI_ub_td_emp_List = zeros(n_mods,nt);
shuffmean_ub_emp = zeros(nt,1);
shuffstd_ub_emp = zeros(nt,1);
shuffsig_ub_emp = zeros(nt,1);



for R = 1:n_mods

    f = flist{R};


    %calculate firing rates according to log-normal distribution
    %obtain timd-dependent mean and standard deviation of the log-rate, which
    %follow a normal distribution

    %change path of load directory to make sure it points to correct data
    tmpstr = split(f.ops.dataname,'/');
    f.ops.dataname = strcat([ops.datadir,tmpstr{end}]);        

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

    %parse up trials and simulate fro model
    %log-normal one: lambda = exp(mu_t + sigma^2/2)
    [dat,lambda_ln,~] = cleantrials(f,trialmask,ops.alignment, ops.tmin,ops.tmax,struct('simtype','lognormal'));
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
    dat = dat(condmask_all,:);
    

    %calculate spiking entropy from the empirical spike count distriution
    %assume max of 10 spikes/bin (200 Hz firing)

    %% the upper limit of mutual information from Jensen's inequality
    maxct = 10;
    H = poiss_entropy(nanmean(lambda_ln)*f.ops.dt,maxct); %should this be from data or model?
    H_glm = poiss_entropy(nanmean(lambda_glm)*f.ops.dt,maxct); %glm form
    H_emp =  poiss_entropy(nanmean(dat)*f.ops.dt,maxct); %emperical form

    %TODO: Use either ops.IsCondScalar or iscell() to decide if you should
    %extract from ops.condmask
       
    n_cond = numel(ops.condmask);
    [ns,nt] = size(lambda_ln);

    %ntrial = numel(f.ops.Dat.trial_idx);

    pvec = zeros(n_cond,1); %probability of stimulus event occuring. use emperical distribution
    p_lam_vec = zeros(n_cond,nt); %mean of log-normal, conditioned on stim. config.
    H_j = zeros(n_cond,nt);
    
    p_lam_glm_vec = zeros(n_cond,nt); %mean of glm, conditioned on stim. config.
    H_j_glm = zeros(n_cond,nt);
    
    %empirical form
    p_lam_emp_vec = zeros(n_cond,nt); %mean of glm, conditioned on stim. config.
    H_j_emp = zeros(n_cond,nt);

    for j = 1:n_cond
        %mask out each condition as categories of the stimulus distribution
        %account for edge clipping
        condtrial = f.ops.Dat.trial_idx(ops.condmask{j});
        condmask = ismember(trial_idx_clip,condtrial);

        pvec(j) = sum(condmask)/ns;

        p_lam_vec(j,:) = nanmean(lambda_ln(condmask,:));
        p_lam_glm_vec(j,:) = nanmean(lambda_ln(condmask,:));
        p_lam_emp_vec(j,:) = nanmean(dat(condmask,:));
        

        H_j(j,:) = poiss_entropy(p_lam_vec(j,:)*f.ops.dt,maxct);
        H_j_glm(j,:) = poiss_entropy(p_lam_glm_vec(j,:)*f.ops.dt,maxct);
        H_j_emp(j,:) = poiss_entropy(p_lam_emp_vec(j,:)*f.ops.dt,maxct);
    end

    Hcond = pvec'*H_j;
    
    MI_lltd = H - Hcond;
    MI_lltd = MI_lltd/log(2);
    MI_ub_tdList(R,:) = MI_lltd;
    
    Hcond_glm = pvec'*H_j_glm;
    MI_lltd_glm = H_glm - Hcond_glm;
    MI_lltd_glm = MI_lltd_glm/log(2);
    MI_ub_td_glm_List(R,:) = MI_lltd_glm;
    
    Hcond_emp = pvec'*H_j_emp;
    MI_lltd_emp = H_emp - Hcond_emp;
    MI_lltd_emp = MI_lltd_emp/log(2);
    MI_ub_td_emp_List(R,:) = MI_lltd_emp;
    
    
    %% create a shuffled version of this to establish a baseline use the upper-bounded
    % point estimate form, since it is fast. just to see
    
     %only do the shuffle on last model, assumed to be the most
     %  sophisticated ones
    if R == n_mods                    
        H = poiss_entropy(nanmean(lambda_ln)*f.ops.dt,maxct);
        H_glm = poiss_entropy(nanmean(lambda_glm)*f.ops.dt,maxct);
        H_emp = poiss_entropy(nanmean(dat)*f.ops.dt,maxct);

        shuff_trials = zeros(nshuff,ns);
        for j = 1:nshuff
            shuff_trials(j,:) = randperm(ns);
        end

        MI_shuff = zeros(nshuff,nt);
        MI_shuff_glm = zeros(nshuff,nt);
        MI_shuff_emp = zeros(nshuff,nt);

        %shufle the trial_idx_clip to scramble the condmask
        for k = 1:nshuff

            trial_idx_clip_shuff = trial_idx_clip(shuff_trials(k,:));
            %ntrial = numel(f.ops.Dat.trial_idx);

            pvec = zeros(n_cond,1); %probability of stimulus event occuring. use emperical distribution
            p_lam_vec = zeros(n_cond,nt); %mean of log-normal, conditioned on stim. config.
            p_lam_glm_vec = zeros(n_cond,nt); %mean of glm, conditioned on stim. config.
            p_lam_emp_vec = zeros(n_cond,nt); %mean of glm, conditioned on stim. config.
            
            H_j = zeros(n_cond,nt);
            H_j_glm = zeros(n_cond,nt);
            H_j_emp = zeros(n_cond,nt);

            for j = 1:n_cond
                %mask out each condition as categories of the stimulus distribution
                %account for edge clipping
                condtrial = f.ops.Dat.trial_idx(ops.condmask{j});
                condmask = ismember(trial_idx_clip_shuff,condtrial);

                pvec(j) = sum(condmask)/ns;

                p_lam_vec(j,:) = nanmean(lambda_ln(condmask,:));
                p_lam_glm_vec(j,:) = nanmean(lambda_glm(condmask,:));
                p_lam_emp_vec(j,:) = nanmean(dat(condmask,:));

                H_j(j,:) = poiss_entropy(p_lam_vec(j,:)*f.ops.dt,maxct);
                H_j_glm(j,:) = poiss_entropy(p_lam_glm_vec(j,:)*f.ops.dt,maxct);
                H_j_emp(j,:) = poiss_entropy(p_lam_emp_vec(j,:)*f.ops.dt,maxct);

            end

            Hcond = pvec'*H_j;
            Hcond_glm = pvec'*H_j_glm;
            Hcond_emp = pvec'*H_j_emp;
            

            MI_shuff(k,:) = H - Hcond;
            MI_shuff_glm(k,:) = H_glm - Hcond_glm;
            MI_shuff_emp(k,:) = H_emp - Hcond_emp;
        end

        shuffmean_ub = nanmean(MI_shuff);
        shuffstd_ub = nanstd(MI_shuff);
        
        shuffmean_ub_glm = nanmean(MI_shuff_glm);
        shuffstd_ub_glm = nanstd(MI_shuff_glm);
        
        shuffmean_ub_emp = nanmean(MI_shuff_emp);
        shuffstd_ub_emp = nanstd(MI_shuff_emp);
        
        %for every time point, calculate the 95% CI
        for l = 1:nt  
            mh = MI_shuff(~isnan(MI_shuff(:,l)),l);
            if ~isempty(mh)
                [h,x] = histcounts(mh,20,'Normalization','cdf');
                g95 = find(h > 0.975);
                shuffsig_ub(l) = x(g95(1)); %first instance
            end
            
            mh = MI_shuff_glm(~isnan(MI_shuff_glm(:,l)),l);
            if ~isempty(mh)
                [h,x] = histcounts(mh,20,'Normalization','cdf');
                g95 = find(h > 0.975);
                shuffsig_ub_glm(l) = x(g95(1)); %first instance
            end
            
            mh = MI_shuff_emp(~isnan(MI_shuff_emp(:,l)),l);
            if ~isempty(mh)
                [h,x] = histcounts(mh,20,'Normalization','cdf');
                g95 = find(h > 0.975);
                shuffsig_ub_emp(l) = x(g95(1)); %first instance
            end
        end
        

    end
       
    %package together all of the glm stuff
    MI_glm_ub = struct();
    MI_glm_ub.MI = MI_ub_td_glm_List;
    MI_glm_ub.shuffmean = shuffmean_ub_glm;
    MI_glm_ub.shuffstd = shuffstd_ub_glm;
    MI_glm_ub.shuffsig = shuffsig_ub_glm;
    
    %package together all of the empirical stuff
    MI_emp = struct();
    MI_emp.MI = MI_ub_td_emp_List;
    MI_emp.shuffmean = shuffmean_ub_emp;
    MI_emp.shuffstd = shuffstd_ub_emp;
    MI_emp.shuffsig = shuffsig_ub_emp;

    %% the sampling approach
    
    shuff_trials = zeros(nshuff,ns);
    for j = 1:nshuff
        shuff_trials(j,:) = randperm(ns);
    end
    

    %basics
    maxct = 10;
    n_cond = numel(ops.condmask);
    [ns,nt] = size(lambda_ln);
    pvec = zeros(n_cond,1);

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
    %check to see if done correctly
    %{
    clf
    hold on
    plot(lambda_ln(1,:))
    plot(nanmean(squeeze(lambda_samp(:,1,:))))
    %}

    %propagate to poisson, and build count distribution
    %counts with rate parameter drawn from above samples of lambda
    y_samps = zeros(nsamp,ns,nt);
    for j = 1:ns
        for k = 1:nt
            y_samps(:,j,k) = poissrnd(lambda_samp(:,j,k)*f.ops.dt);
        end
    end
    %check
    %plot(nanmean(squeeze(y_samps(:,1,:)/0.05)))

    %create the conditional distributions over spike counts
    y_counts = zeros(maxct,ns,nt);
    for j = 1:ns
        for k = 1:nt
            y_counts(:,j,k) = histcounts(y_samps(:,j,k),0:maxct);
        end
    end
    %check
    %{
    spike_vector = 0:maxct-1;
    test = squeeze(y_counts(:,1,:));
    test = test./sum(test,1);
    meanct = spike_vector*test;
    plot(meanct/0.05)
    %}


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
        
        %check
        %{
        spike_vector = 0:maxct-1;
        test = p_j;
        meanct = spike_vector*test;
        plot(meanct/0.05,'linewidth',2)
        %}

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

    
    
    %check
     %{
    spike_vector = 0:maxct-1;
    test = p_marg;
    meanct = spike_vector*test;
    
    plot(meanct/0.05,'linewidth',2)
     %}

    
    MI_tdList(R,:) = MI_td;

    %% figure plotting

    if makeFig
        plot(tmesh,MI_td,'color',colormat(R,:),'linewidth',2);
        %plot(tmesh,MI_lltd,'--','color',colormat(R,:),'linewidth',2);  
        xlabel(strcat(['time to ',ops.alignment,' (s)']));
        ylabel('MI')
        title(strcat(['MI (bits), Neuron ',num2str(file_idx),' stim: ',ops.stimname]));
        set(gca,'fontsize',15);
    end
    
    
end

if makeFig
    
    %originally had all models iwth their own shuffle mean. ignore
    %{
    for k = n_mods %just plot final one
        %shuffmean = nanmean(MI_ub_shuffList(R,:,:));
        %shuffsig = 2*nanstd(MI_shuff(R,:,:));
        shuffmean_samp = nanmean(MI_shuffList(k,:,:));
        shuffsig_samp = 2*nanstd(MI_shuffList(k,:,:));
        
        %shadedErrorBar(tmesh,shuffmean,shuffsig,'lineprops',{'color','blue'})
        shadedErrorBar(tmesh,shuffmean_samp,shuffsig_samp) 
    end
    %}
    shadedErrorBar(tmesh,shuffmean,2*shuffstd); 
    
    vline(0,'k')
    
    %decide if you want to show upper bound or not
    %llist = cell(2*n_mods,1);
    llist = cell(n_mods+1,1);
    ind = 0;
    for k = 1:n_mods
        llist{ind+1} = strcat([ops.modelnames{k}]);
        %llist{ind+2} = strcat(['point-est:',ops.modelnames{k}]);
        %ind = ind+2;
        ind = ind + 1;
    end
    llist{ind+1} = 'shuffle';
    legend(llist,'location','best','Interpreter','None')

    %}
    
    
    %legend({ops.modelnames{1:n_mods},'shuffle'},'location','best')
end










