function [rate_dat,lambda,trial_idx_final,lambda_var,data_idx,mod_idx,XKnew] = ...
        cleantrials(f,trialmask,alignment, tmin,tmax,ops)
%parse up continuous data into trials and trim edges correctly to compare data and a model fit
%
%Inputs:
    %   f: loaded fit structure
    %   trialmask: (ns x 1) boolean for which trials to simulate
    %   alignment: event around which to center trials: 'start','leavecpoke','choice'
    %   tmin: min time relative to start a trial
    %   tmax: max time after start for a trial
    %   ops: extra options structure with following info
    %       datadir: where raw data files live
    %       simtype: "glm" or "lognormal" to control if simulating from
    %       pointestimate-style GLM, or log-normal distribution GAM
    %       stim: if supplied, a stim struct to replace default stim
    %       
    
    %
    %Outputs, for ns trial, nt time points
    %   rate_dat: (ns x nt) matrix of raw rate data
    %   lambda: (ns x nt) matrix of model-predicted rate
    %   trial_idx_final: clipped version of trial ids, which may be shorter
    %           if edge-clipping removed an impartial trial
    %   lambda_var: (ns x nt) variance of model-predicted rate
    %           ("lognormal" option in ops.simtype only)
    %   data_idx: (ns x nt) indexing matrix for how each lambda entry
    %   relates to raw data
    %   mod_idx: similar for data
    %   XKnew: the contributions from each kernel

    
%can be artifacts of fitting on cluster that stay in f.ops.dataname
%make sure f.ops.fname is pointing to correct dir for these sims. 
if isfield(ops,'datadir')
    tmpstr = split(f.ops.dataname,'/');
    f.ops.dataname = strcat([ops.datadir,tmpstr{end}]);        
end

%determine which params to use, make stim and dat structures
[~,kuse] = min(f.NLL.val); %which kfold to use

[Dat,stim] = parse_stimuli_wrapper(f.ops);

%overwrite default stimuli with supplied one
if isfield(ops,'stim')
    stim = ops.stim;
end

if isfield(ops,'trial_idx')
    Dat.trial_idx = ops.trial_idx; %way to ovewrite trials
end

%TEMP REMOVE: updates ops sizing if a special form of basis has been used
%if strcmp(f.ops.basistype,'ofc_adaptive_plus')
%    f.ops.dimst(:) = 9;   
%end

f.ops.Dat = Dat;
[basefun,basefunst] = chooseBasis_ofc(f.ops);
stim.basefuns = {basefun,basefunst};

%TEMP: updates ops sizing if a special form of basis has been used
%if strcmp(f.ops.basistype,'ofc_adaptive_plus')
%    nbnew = size(stim.basefuns{2}{1},1);
%    f.ops.dimst(:) = nbnew;   
%end

rate_datall = smooth(Dat.y/f.ops.dt); %correct smoothing

%trial numbers of trialmask. references original data, which includes all
%data (e.g., opt-out trials included).
%the model data, on the other hand, will have had that extra data removed
%and has a different reference frame for indexing
trial_idx = Dat.trial_idx(trialmask);



%define minimum and max of area
Tmin = tmin/f.ops.dt; %2s before align_idx
Tmax = tmax/f.ops.dt; %4s after enter cpoke


%simulate data. make better switch case. But simulate all trials, then filter
if ~isfield(ops,'simtype')
    [~,lambda,timingGLM_idx,XK] = simulateGLM(f.wML{kuse},stim,f.ops,trialmask);
    lambda_var = [];
else   
    if strcmp(ops.simtype,'none') %don't simulate, just grab dat
        lambda =  rate_datall;
        lambda_var = [];
        %essentially duplicate the Dat info
        timingGLM_idx = struct();
        timingGLM_idx.entercpoke_idx = Dat.entercpoke_idx(trial_idx);
        timingGLM_idx.leavecpoke_idx = Dat.leavecpoke_idx(trial_idx);
        timingGLM_idx.choice_idx =  Dat.choice_idx(trial_idx);        
  
    elseif ~strcmp(ops.simtype,'lognormal')
        [~,lambda,timingGLM_idx,XK] = simulateGLM(f.wML{kuse},stim,f.ops,trialmask);
        lambda_var = [];
    else %simulate via log normal
        [~,lambda,timingGLM_idx,lambda_var] = simulateGLM_logNormal(f.wML{kuse},f.wML_cov{kuse},stim,f.ops,trialmask);
        XK = [];
    end
end



% align to what event?
switch alignment
    case 'start'
        align_idx_dat = Dat.entercpoke_idx(trial_idx)';
        align_idx_mod = timingGLM_idx.entercpoke_idx';
    case 'leavecpoke'
        align_idx_dat = Dat.leavecpoke_idx(trial_idx)';
        align_idx_mod = timingGLM_idx.leavecpoke_idx';
    case 'choice'
        align_idx_dat = Dat.choice_idx(trial_idx)';
        align_idx_mod = timingGLM_idx.choice_idx';
    case 'flash'
        %helps parse what to keep for data
        [~,keepmask_dat] = parse_pseudotrials(f.ops,trialmask); 
        align_idx_dat = find(stim.TDstim{strcmp('LFlash',stim.stimlegend)}(keepmask_dat)+...
            stim.TDstim{strcmp('RFlash',stim.stimlegend)}(keepmask_dat));
        align_idx_mod = unique(sort([timingGLM_idx.lflash_idx',timingGLM_idx.rflash_idx']));
    case 'click'
        %helps parse what to keep for data
        [~,keepmask_dat] = parse_pseudotrials(f.ops,trialmask); 
         align_idx_dat = find(stim.TDstim{strcmp('LClick',stim.stimlegend)}+...
            stim.TDstim{strcmp('RClick',stim.stimlegend)});
        align_idx_mod = unique(sort([timingGLM_idx.lclick_idx',timingGLM_idx.rclick_idx']));
end

%align_idx_dat = align_idx_dat(trial_idx);
%align_idx_mod = align_idx_mod(trialmask);

% create data_idx mask to grab data from any model simulation
ns = numel(align_idx_dat);
nt = numel(Tmin:Tmax);
data_idx = zeros(ns,nt); %mask for trial-structured, aligned raw data
mod_idx = zeros(ns,nt);
for j = 1:ns
    data_idx(j,:) = align_idx_dat(j)+Tmin:1:align_idx_dat(j)+Tmax;
    mod_idx(j,:) = align_idx_mod(j) + Tmin:align_idx_mod(j) + Tmax;
end

%check edge cases, but omit values instead of padding with last values
trial_idx_final = trial_idx;
%check the model data overlap with edge cases. if need to omit, do to both
%data and model indices to maintain consistent sizes
%end of trial
for j = numel(align_idx_mod):-1:1
    edge_num = (align_idx_mod(j) + Tmax) - numel(lambda);
    edge_num_dat = (align_idx_dat(j) + Tmax) - numel(rate_datall);
    if (edge_num > 0) || (edge_num_dat > 0)
        mod_idx = mod_idx(1:j-1,:);
        data_idx = data_idx(1:j-1,:);
        trial_idx_final = trial_idx_final(1:end-1);
    else
        break;
    end
end

%beginning of trial
for j = 1:numel(align_idx_mod)
    edge_num = (align_idx_mod(j) + Tmin);
    edge_num_dat = (align_idx_dat(j)  + Tmin);
    if (edge_num < 1) || (edge_num_dat < 1)
        mod_idx = mod_idx(j+1:end,:);
        data_idx = data_idx(j+1:end,:);
        trial_idx_final = trial_idx_final(2:end);
    else
        break;
    end

end

%%aligned, smoothed data firing rate.
%rate_datall = smooth(Dat.y/f.ops.dt); %correct smoothing
rate_dat = rate_datall(data_idx); 

%smooth and grab balanced trials
lambda= smooth(lambda);
lambda = lambda(mod_idx);

%smooth and grab important parts of XK
XKnew = zeros(size(lambda,1),size(lambda,2),size(XK,2));
for j = 1:size(XKnew,3)
    XKtemp = smooth(XK(:,j));
    XKnew(:,:,j) = XKtemp(mod_idx);
end


%handle variance if using log normal
if isfield(ops,'simtype')
    if strcmp(ops.simtype,'lognormal')
        lambda_var = lambda_var(mod_idx);
    end
end


