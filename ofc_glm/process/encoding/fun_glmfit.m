function [NLL,evidence,ops,wML,wML_cov,hyperparam,paramav] = fun_glmfit(ops,varargin)
%core worker function to fit GLMS
%updated version
%
%inputs:
%   ops: data structure containing all important aspects of training
%ouputs:
%   NLL: negative log likelihood of optimal solution
%   evidence: value of evidence, useful for hyperparameter fits

%% format input data
fname = ops.fname; %eventual save name

%set verobisty level to max
if ~isfield(ops,'dispL')
    ops.dispL = 5;
end

% load data (if not loaded already) and parse stimuli into stimulus struct
dispL('parse stimuli',3,ops.dispL);
if isempty(varargin) %stim not supplied
    [Dat,stim] = parse_stimuli_wrapper(ops);
    ops.Dat = Dat;
else
    [Dat,~] = parse_stimuli_wrapper(ops);
    ops.Dat = Dat;
    stim = varargin{1};
end
ops.nkern = sum(stim.TDcausalvec,'all'); %captures changes due to omit kernels
%ops.nkern = Dat.nkern;


%create basis functions
dispL('basis functions',3,ops.dispL)
[basefun,basefunst] = chooseBasis_ofc(ops);
stim.basefuns = {basefun,basefunst};

%TEMP REMOVE: updates ops sizing if a special form of basis has been used
%if strcmp(ops.basistype,'ofc_adaptive_plus')
%    nbnew = size(stim.basefuns{2}{1},1);
%    ops.dimst(:) = nbnew;   
%end

%separate out y for ease. TODO: omit for memory?
y = ops.Dat.y; %spike times. size (nt,ns). nt=time points, ns = number trials, even if pseudotrial
nt = size(y,1);
ops.nt = nt; %TODO: is this needed?

%decide on how to parition trials into stratified cross-validated parts
kfold = ops.kfold; %number of cross validation parts
kfold_inds = cvpartition_wrapper(ops);
ops.kfold_inds = kfold_inds;
kmax = max(kfold-1,1); %if kfold==1, it's one whole data set. no test set

%saved data structures after final hyperparam is chosen.

%NLL_test = zeros(1,kmax); %negative log-likelihood for test set
%evidence_test = zeros(1,kmax); %log evidence for test set
%NLL_val = zeros(1,kmax); %negative log-likelihood for validation set
%evidence_val = zeros(1,kmax); %log evidence for validation set
%NLL_train = zeros(1,kmax); %negative log-likelihood for training set
%evidence_train = zeros(1,kmax); %log evidence for training set

NLL = struct();
evidence = struct();
%NLL=struct('test',zeros(1,kmax), 'train',zeros(1,kmax),'val',zeros(1,kmax));
%evidence=struct('test',zeros(1,kmax), 'train',zeros(1,kmax),'val',zeros(1,kmax));

wML = cell(1,kmax); %optimal params, per kfold
wML_cov = cell(1,kmax); %covariance of parameters, per kfold. based on test set


%base penalty terms, not weighted by number of samples
%gamma_base = ops.gamma;
%beta_base = ops.beta;
hyperparam = struct(); %struct to hold information about 
                            %hyperparam grid search

%num controls.
if ~isempty(stim.TDcausalvec)
    TDcaus_inds = find(stim.TDcausalvec(1,:)); %indices, causal filtering
    TDacaus_inds = find(stim.TDcausalvec(2,:)); %indices, acausal filtering
else
    TDcaus_inds = []; 
    TDacaus_inds = [];
end
nx = numel(stim.xmat); %# time-independent parameters, including bias
dim_stim = sum(ops.dimst(TDcaus_inds))+sum(ops.dimst(TDacaus_inds)); %total # basis functions for stimuli
nc = sum(ops.dim) + dim_stim + nx; %total number of params

%the number of time points for each fold
ops.nsnt = zeros(kfold,1);

%decide on algorithm type. L1 can be weird with trust-point, but seems
%stable with quasi-newton
if strcmp(ops.regularizer,'L2')
    algtype = 'trust-region'; %fast, but can be weird for L1
elseif strcmp(ops.regularizer,'L1')
    algtype = 'quasi-newton'; %slower, but not weird for L1
elseif strcmp(ops.regularizer,'L1+L2')
    algtype = 'quasi-newton'; %slower, but not weird for L1
end

%params for grid search over hyperparameter values
%default grid search setup. Also decide on alg type. if L1 included, use
%quasi newton. otherwise, try trust point (faster)
if ops.useHyperParamGrid
    dispL('Performing hyperparam search',3,ops.dispL)
    n_reglevels = 4; %number levels to refine search
    n_regmesh = 5; %number of values on grid at each level to search
    
    %grid to do original sampling. save output for analysis  
    %regularizer_initvec = [1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-2,1e-1,1,10,100];
    %regularizer_initvec = [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1,10,100];
    
    regularizer_initvec = [1e-5,1e-4,1e-3,1e-2,1e-1,1,10];
    
    %for first grid level
    %start with default values either given by user; or a
    %default, log-linear set to cast a wide net
    %TODO: eventually generalize to any number and type of penalties
    if strcmp(ops.regularizer,'L2')
            gamma_vec = regularizer_initvec;
            beta_vec = 0;
            npenalty = 1;
    elseif strcmp(ops.regularizer,'L1')
            beta_vec = regularizer_initvec;
            gamma_vec = 0;
            npenalty = 1;           
    elseif strcmp(ops.regularizer,'L1+L2') %use both
            gamma_vec = regularizer_initvec;
            beta_vec = regularizer_initvec;
            npenalty = 2;
    end

else
    dispL('Bypassing hyperparam search',3,ops.dispL)
    n_reglevels = 1;
    n_regmesh = 1;
    if strcmp(ops.regularizer,'L2')
            gamma_vec = ops.gamma;
            beta_vec = 0;
            npenalty = 1;
    elseif strcmp(ops.regularizer,'L1')
            beta_vec = ops.beta;
            gamma_vec = 0;
            npenalty = 1;
    elseif strcmp(ops.regularizer,'L1+L2') %use both
            gamma_vec = ops.gamma;
            beta_vec =  ops.beta;
            npenalty = 2;
    end
        
end

%loop over levels of regularizer testing
%loop over the grid values
for gridlevel = 1:n_reglevels       
    if ops.useHyperParamGrid
        dispL(strcat(['hyperparam search gridlevel:',num2str(gridlevel)]),3,ops.dispL)
    end

    if gridlevel > 1
        %refined testing values for gamma and beta         
        if strcmp(ops.regularizer,'L2')
            gamma_vec = linspace(gamma_1,gamma_2,n_regmesh);
            beta_vec = 0;
        elseif strcmp(ops.regularizer,'L1')
            beta_vec = linspace(beta_1,beta_2,n_regmesh);
            gamma_vec = 0;
        elseif strcmp(ops.regularizer,'L1+L2') %use both
            beta_vec = linspace(beta_1,beta_2,n_regmesh);
            gamma_vec = linspace(gamma_1,gamma_2,n_regmesh);
        end
    end
    
    %results of grid search go here
    evidence_mat_test = zeros(numel(beta_vec),numel(gamma_vec),kmax); %for test set
    evidence_mat_val = zeros(numel(beta_vec),numel(gamma_vec),kmax); %for validation set
    evidence_mat_train = zeros(numel(beta_vec),numel(gamma_vec),kmax); %for training set
    
    nll_mat_test = zeros(numel(beta_vec),numel(gamma_vec),kmax); %for test set
    nll_mat_val = zeros(numel(beta_vec),numel(gamma_vec),kmax); %for validation set
    nll_mat_train = zeros(numel(beta_vec),numel(gamma_vec),kmax); %for training set
    
    %pillow type nll. lacks a dc term
    dc_mat_test = zeros(numel(beta_vec),numel(gamma_vec),kmax); %for test set
    dc_mat_val = zeros(numel(beta_vec),numel(gamma_vec),kmax); %for validation set
    dc_mat_train = zeros(numel(beta_vec),numel(gamma_vec),kmax); %for training set
    
    %regularizer values
    penalty_mat = zeros(numel(beta_vec),numel(gamma_vec),kmax,npenalty);
    
    %values for use in hyperparameter optimization
    evidence_mat_hyp = zeros(numel(beta_vec),numel(gamma_vec));
    nll_mat_hyp = zeros(numel(beta_vec),numel(gamma_vec));
    
    %which kfold was chosen;
    kusevec = zeros(numel(beta_vec),numel(gamma_vec)); 
    
    %use to later propagate final param values for each k-fold
    wML_mat_all = zeros(numel(beta_vec),numel(gamma_vec),kmax,nc);
        
    %start with "Piggyback" off to use a new initial seed for params. turn
    %on after first optimization
    keepPiggyback = false;
    
    %use numel() instead of n_regmesh to short circuit loop if not
    %using a given penalty term
    for beta_ind = 1:numel(beta_vec)

        for gamma_ind = 1:numel(gamma_vec)

            if ops.useHyperParamGrid
                dispL(strcat(['meshpoint [',...
                    num2str(beta_ind),', ',num2str(gamma_ind),']']),3,ops.dispL)
            end
       
            ops.test_ids = kfold_inds == kfold; %the overall held-out test ids
            
            for kfsamp = 1:kmax

                if kfold > 1
                    ops.val_ids = kfold_inds==kfsamp;
                    %exclude last kfold. test data
                    ops.train_ids = (kfold_inds ~=kfsamp) & (kfold_inds ~= kfold); 
                else %no paritioning. probably not advisable, but good test condition
                    ops.val_ids = kfold_inds==kfsamp;
                    ops.train_ids = ops.test_ids;
                end   


                %find number of data points to normalize LL values
                train_id = ops.Dat.trial_idx(ops.train_ids);
                val_id = ops.Dat.trial_idx(ops.val_ids);
                test_id = ops.Dat.trial_idx(ops.test_ids);              
                nsnt_train = sum(range([ops.Dat.start_idx(train_id);ops.Dat.end_idx(train_id)])+1);
                nsnt_val = sum(range([ops.Dat.start_idx(val_id);ops.Dat.end_idx(val_id)])+1);
                nsnt_test = sum(range([ops.Dat.start_idx(test_id);ops.Dat.end_idx(test_id)])+1);
                
                %set number of timepoints
                ops.nsnt(kfsamp) = nsnt_val;
                ops.nsnt(kfold) = nsnt_test;
                              
                ops.beta = beta_vec(beta_ind);
                ops.gamma = gamma_vec(gamma_ind);


                %% fit the GLM and save
                

                if isfield(ops,'prevname') %daisy-chain a prevous trial's results 
                    %TODO: REVISE THIS CODE. NO LONGER USED THAT OFTEN
                    disp('loading previous file for optimization params');
                    a = load(ops.prevname);
                    if isfield(a,'wML')
                        wML_prev = a.wML;
                    else %a temporary sav file
                        disp('opening from a temporary save file')
                        wML_prev = a.x;
                    end

                    %deal out previous params.        
                    ci = -rand(1,nc);
                    if numel(wML_prev) < nc
                        ci(1:numel(wML_prev)-1) = wML_prev(1:end-1);
                    end
                    ci(end) = wML_prev(end);
                elseif keepPiggyback %use a previous run in the kfold + hyperparam loop. its a convex optimization
                    ci = wML_k;
                else
                    dispL('using new initial point for optimiation',5,ops.dispL);
                    ci = -rand(1,nc); %choose - vals for a defined initial LL point
                end

                %do an initial check
                %if ~ops.useHyperParamGrid
                    dispL('initial normed LL: training data',4,ops.dispL);
                    %dispL(1/(sf)*nll_poiss(y,ci,stim,ops,ops.train_ids),1,ops.dispL);
                    dispL(1/nsnt_train*nll_poiss(y,ci,stim,ops,ops.train_ids),4,ops.dispL);
                    
                    dispL('initial normed LL: validation data',4,ops.dispL);
                    %dispL(1/(sf_test)*nll_poiss(y,ci,stim,ops,ops.train_ids),1,ops.dispL);
                    dispL(1/nsnt_val*nll_poiss(y,ci,stim,ops,ops.val_ids),4,ops.dispL);

                    dispL('initial normed LL: test data',4,ops.dispL);
                    %dispL(1/(sf_test)*nll_poiss(y,ci,stim,ops,ops.train_ids),1,ops.dispL);
                    dispL(1/nsnt_test*nll_poiss(y,ci,stim,ops,ops.test_ids),4,ops.dispL);
                %end

                % fit parameters
                if ~isfield(ops,'funiter')
                    maxiter = 2000;
                else
                    maxiter = ops.funiter;
                end

                %decide on displaying output
                if ops.dispL > 4
                   useDisplay = 'iter';
                else 
                    useDisplay = 'off';
                end

                %optimize
                if ~ops.useHyperParamGrid
                    dispL('optimizing',1,ops.dispL)
                end
                if ~ops.isConstrained
                    myoutput_fcn = @(x,optimvalues,state) myoutput(x,optimvalues,state,ops);
                    options = optimoptions('fminunc','GradObj','on',...
                    'Hessian','on','algorithm',algtype,'display',useDisplay,...
                    'maxfunevals',8000,'maxIterations',maxiter,'functionTolerance',1e-6,...
                    'OutputFcn',myoutput_fcn);

                    func = @(x) nll_poiss(y,x,stim,ops,ops.train_ids);
                    [wML_k, ~, ~, ~,~, ~] = fminunc(func, ci, options);

                else
                    %try restricting the history filter to be strictly negative   
                    lb = -inf*ones(nc,1);
                    ub = inf*ones(nc,1); ub(1:ops.dim) = 0; 
                    myoutput_fcn = @(x,optimvalues,state) myoutput(x,optimvalues,state,ops);
                    func = @(x) nll_poiss(y,x,stim,ops,ops.train_ids);
                     options = optimoptions('fmincon','GradObj','on',...
                    'Hessian','on','algorithm','trust-region-reflective',useDisplay,'iter',...
                    'maxfunevals',8000,'maxIterations',maxiter,'functionTolerance',1e-6,...
                    'OutputFcn',myoutput_fcn);
                    %for regular constraint to keep negative
                    [wML_k, ~, ~, ~,~, ~] = fmincon(func,ci,[],[],[],[],lb,ub,[],options);

                end
                
                %turn on piggyback after first run
                keepPiggyback = true;
                
                % output for a validation and training data set     
                [NLL_ktrain,~,~,evidence_ktrain,penalty_k,dc_train] = nll_poiss(y,wML_k,stim,ops,ops.train_ids);   
                %validation set nLL  and evidence                     
                [NLL_kval,~,~,evidence_kval,~,dc_val] = nll_poiss(y,wML_k,stim,ops,ops.val_ids);
                %test set nLL  and evidence                     
                [NLL_ktest,~,~,evidence_ktest,~,dc_test] = nll_poiss(y,wML_k,stim,ops,ops.test_ids);
                
                %if ~ops.useHyperParamGrid
                    dispL(fname,4,ops.dispL)
                    
                    dispL('normed training set nLL  and evidence',4,ops.dispL);          
                    %dispL([1/sf*NLL_k, 1/sf*evidence_k],1,ops.dispL);
                    dispL(1/nsnt_train*[NLL_ktrain, evidence_ktrain],4,ops.dispL);

                    dispL('normed validation set nLL  and evidence.',4,ops.dispL);
                    dispL(1/nsnt_val*[NLL_kval, evidence_kval],4,ops.dispL);
                    
                    dispL('normed testing set nLL and evidence.',4,ops.dispL);
                    dispL(1/nsnt_test*[NLL_ktest, evidence_ktest],4,ops.dispL);
                %end                
                
                %store train,val, test 
                nll_mat_test(beta_ind,gamma_ind,kfsamp) = NLL_ktest;
                nll_mat_val(beta_ind,gamma_ind,kfsamp) = NLL_kval;
                nll_mat_train(beta_ind,gamma_ind,kfsamp) = NLL_ktrain;
                
                %TODO: EVENTUALLY REMOVE
                dc_mat_test(beta_ind,gamma_ind,kfsamp) = dc_test;
                dc_mat_val(beta_ind,gamma_ind,kfsamp) = dc_val;
                dc_mat_train(beta_ind,gamma_ind,kfsamp) = dc_train;
                
                evidence_mat_test(beta_ind,gamma_ind,kfsamp) = evidence_ktest;
                evidence_mat_val(beta_ind,gamma_ind,kfsamp) = evidence_kval;
                evidence_mat_train(beta_ind,gamma_ind,kfsamp) = evidence_ktrain;               
                
                %keep all param values for hyperparam,kfold combo
                wML_mat_all(beta_ind,gamma_ind,kfsamp,:) = wML_k;
                
                %penalty term
                penalty_mat(beta_ind,gamma_ind,kfsamp,:) = penalty_k;             

                %save a temp version
                %save(strcat([fname,'_ktemp.mat']),'wML','stim*','ops','nll_mat*','evidence_mat*','wML_cov','hyperparam')

            end
    
            
            %best kfold test set nLL and evidence values for current hyperparam fit
            [~,k_ind_best] = min(nll_mat_val(beta_ind,gamma_ind,:));
            kusevec(beta_ind,gamma_ind) = k_ind_best;       
            
            %values for use in hyperparameter optimization
            %%use log evidence
            %evidence_mat_hyp(beta_ind,gamma_ind) = evidence_mat_test(beta_ind,gamma_ind,k_ind_best);       
            %%log-evidence of training data (might be correct?)
            %evidence_mat_hyp(beta_ind,gamma_ind) = evidence_mat_train(beta_ind,gamma_ind,k_ind_best);
            
            %hacky way to use test likelihood
            evidence_mat_hyp(beta_ind,gamma_ind) = ...
                -(nll_mat_test(beta_ind,gamma_ind,k_ind_best)-sum(squeeze(penalty_mat(beta_ind,gamma_ind,k_ind_best,:))));
            
            
            
            nll_mat_hyp(beta_ind,gamma_ind) = nll_mat_test(beta_ind,gamma_ind,k_ind_best);
            
            %alternatively, average together parameters and calculate test
            %set
            %{
            wML_temp = zeros(kmax,nc);
            for kfsamp = 1:kmax
                wML_temp(kfsamp,:) = wML{kfsamp};
            end
            wML_av = mean(wML_temp);
            %test set nLL  and evidence                     
            [NLL_ktest,~,~,evidence_ktest] = nll_poiss(y,wML_av,stim,ops,ops.test_ids);
            nll_mat(beta_ind,gamma_ind) = NLL_ktest;
            evidence_mat(beta_ind,gamma_ind) = evidence_ktest;
            %wML_kbest = wML{k_ind_best};
            wML_mat(beta_ind,gamma_ind,:) = wML_av;
            %}

        end
    end
        
    %decide on beta_min/max and gamma min/max for next level of grid
    %search

    %sometimes evidence goes to infinity for really strong penalties. omit.
    evidence_mat_hyp(isinf(evidence_mat_hyp))=nan;

    ev_max = max(evidence_mat_hyp,[],'all');
    [emax_row,emax_col] = find(evidence_mat_hyp==ev_max);
    beta_1 = beta_vec(emax_row);
    gamma_1 = gamma_vec(emax_col);
    [nbeta,ngamma] = size(evidence_mat_hyp); %size of current search grid
    
    %decide if a second value should be found to establish bounds
    %i.e., performign optimization of beta and gamma, or just one
    if numel(beta_vec)==1
        beta_2=beta_1;
    else
        %choose end of range for next iteration as either value to left or
        %right in the grid. 
        %TODO: this can be much cleaner. maybe add during call to ev_max
        if emax_row==nbeta
            beta_2 = beta_vec(end-1);
        elseif emax_row ==1
            beta_2 = beta_vec(2);
        else
            [~,maxind2] = max([evidence_mat_hyp(emax_row-1,emax_col),...
                                evidence_mat_hyp(emax_row+1,emax_col)]);
            if maxind2==1                
                beta_2 = beta_vec(emax_row-1);
            else
                beta_2 = beta_vec(emax_row+1);
            end
        end
        
    end
    if numel(gamma_vec)==1
        gamma_2 = gamma_1; %which should = 0;
    else
        if emax_col==ngamma
            gamma_2 = gamma_vec(end-1);
        elseif emax_col ==1
            gamma_2 = gamma_vec(2);
        else
            [~,maxind2] = max([evidence_mat_hyp(emax_row,emax_col-1),...
                                evidence_mat_hyp(emax_row,emax_col+1)]);
            if maxind2==1                
                gamma_2 = gamma_vec(emax_col-1);
            else
                gamma_2 = gamma_vec(emax_col+1);
            end
        end
        
    end

    %best parameter
    bestgamma = gamma_1;
    bestbeta = beta_1;
    
    %set best regularizers for final output of that kfold  
    ops.beta = bestbeta;
    ops.gamma = bestgamma;
    
    
    %set the values of output data for best grid candidate.
    wML_temp = zeros(kmax,nc);
    for kfsamp = 1:kmax      
        wML{kfsamp} = squeeze(wML_mat_all(emax_row,emax_col,kfsamp,:))';
        
        [~,~,H_k,~] = nll_poiss(y,wML{kfsamp},stim,ops,ops.test_ids);
        wML_temp(kfsamp,:) = wML{kfsamp};
        wML_cov{kfsamp} = inv(H_k);    
        
        %get values of best beta,gamma parameter. what most people will
        %care about.
        NLL.test(kfsamp) = nll_mat_test(emax_row,emax_col,kfsamp);
        NLL.val(kfsamp) = nll_mat_val(emax_row,emax_col,kfsamp);
        NLL.train(kfsamp) = nll_mat_train(emax_row,emax_col,kfsamp);
        
        evidence.test(kfsamp) = evidence_mat_test(emax_row,emax_col,kfsamp);
        evidence.val(kfsamp) = evidence_mat_val(emax_row,emax_col,kfsamp);
        evidence.train(kfsamp) = evidence_mat_train(emax_row,emax_col,kfsamp);
    end

    %create the average parameters and make a sub-structure to record its
    %data
    wML_av = mean(wML_temp);
    [nll_av,~,H_av,evidence_av] = nll_poiss(y,wML_av,stim,ops,ops.test_ids);
    
    paramav = struct();
    paramav.wML = wML_av;
    paramav.NLL = nll_av;
    paramav.evidence = evidence_av;
    paramav.cov = inv(H_av);
    
    
    %do some checks after first pass to see if search should be altered
    if gridlevel==1
        %save original sampling to main file for analysis
        hyperparam.beta_vec = beta_vec;
        hyperparam.gamma_vec = gamma_vec;
        
        hyperparam.nll_test = nll_mat_test;
        hyperparam.evidence_test = evidence_mat_test;
        hyperparam.nll_train = nll_mat_train;
        hyperparam.evidence_train = evidence_mat_train;
        hyperparam.nll_val = nll_mat_val;
        hyperparam.evidence_val = evidence_mat_val;
        hyperparam.penalty = penalty_mat;
        hyperparam.wml = wML_mat_all;
        
        hyperparam.nll_test_dc = dc_mat_test;
        hyperparam.nll_train_dc = dc_mat_train;
        hyperparam.nll_val_dc = dc_mat_val;
        
       
        hyperparam.kuse = kusevec;

        %use heuristic to decide if search will continue.
        %Case to check: if grid search was employed, then gridlevel 1 
        %has more than one entry. if lowest entry was chosen,
        %search will be bogus. set to 1e-8 and move on
        
        %TODO: Once general penalties are coded, clean this up
        switch ops.regularizer
            case 'L1'         
                if numel(beta_vec) > 1 && emax_row==1
                    dispL('hyperparameter search seems fruitless. stopping',3,ops.dispL)
                    break
                end
            case 'L2'
                if numel(gamma_vec) > 1 && emax_col==1
                    dispL('hyperparameter search seems fruitless. stopping',3,ops.dispL)
                    break
                end
            case 'L1+L2'
                if numel(gamma_vec) > 1 && emax_col==1 && numel(beta_vec) > 1 && emax_row==1
                    dispL('hyperparameter search seems fruitless. stopping',3,ops.dispL)
                    break
                end
            otherwise %TODO: need to handle other regularizer setups              
        end
    end   
end
    

      

 % save   
if isfile(fname)
    %remove tmp files
     %delete(strcat([fname,'.tmp']));
     %delete(strcat([fname,'_ktemp.mat']));
     disp('save file already exists. dont overwrite previous data unless you mean it! adding a suffix')
     fname = strcat([fname,'(1).mat']);     
     save(fname,'wML','paramav','stim*','ops','NLL','evidence','wML_cov','hyperparam')

else
    %remove tmp files
    %delete(strcat([fname,'.tmp']));
    %delete(strcat([fname,'_ktemp.mat']));
    save(strcat([fname,'.mat']),'wML','paramav','stim*','ops','NLL','evidence','wML_cov','hyperparam')
end



