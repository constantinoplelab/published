function [stim_new] = remove_covariate(stim,covname)
%removes a covariate from an existing stimulus struct
%inputs:
%   stim: old stim structure for testing/training data
%   covname: name of covariate/kernel to be removed
%
%
%
%outputs: stim_new. new, updated stim struct

stim_new = stim; %initialize

%important sizing information
if isempty(stim.TDcausalvec)
    ntdc = 0;
    ntdac = 0;
else
    ntdc = sum(stim.TDcausalvec(1,:)); %number time-dep. causal stim
    ntdac = sum(stim.TDcausalvec(2,:)); %number time-dep acausal stim
end

nc = ntdc+ntdac;

[nx,~] = size(stim.xmat);
nTDstim = numel(stim.TDstim);

%extract information about location and type of stimuous
%stimind: index in stimlegend
%causal_ind: index in its respective causalvec;

[stimind,causal_ind,covops] = stimulusInfo(stim,covname);
isCausal = covops.isCausal; %causal or acausal stimulus
isConst = covops.isConst; %kernel or constant
keepStim = covops.keepStim; %should stimulus be removed if this kernel is omitted


%remove covariates from either xmat, or TDstim. update stimlegend or
%xorder
if isConst %i.e., a constant bias param
    
    keepinds = find(1:nx~=stimind);
    xmat = stim.xmat(keepinds,:);
    stim_new.xmat = xmat;
    stim_new.Xorder = stim.Xorder(keepinds);

else %a kernel must be removed
    keepinds_all = 1:nc ~=stimind;
    stim_new.stimlegend = stim.stimlegend(keepinds_all);    
        
    if keepStim && isCausal %only set causalvec entry to 0 and move on with life
        stim_new.TDcausalvec(1,causal_ind) = 0;
    elseif keepStim && ~isCausal
        stim_new.TDcausalvec(2,causal_ind) = 0;
    else %remove the stimulus
        keepinds = 1:nTDstim ~=causal_ind;
        stim_new.TDstim = stim.TDstim(keepinds);
        stim_new.TDcausalvec = stim.TDcausalvec(:,keepinds);
    end 
         
end


    
