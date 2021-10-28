function [stimind,causal_ind,covops] = stimulusInfo(stim,covname)
 %grabs index and information about a covariate in stim
 %useful for removing, merging, etc.
 
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

%nSEstim = size(stim.SEstim,2);

%find covariate. first check kernel space
stimind = 0;
isConst = false;
for j = 1:nx
    if strcmp(covname,stim.Xorder{j})
        stimind = j;
        isConst = true;
        break;
    end
end
%then check constant space
for j = 1:nc
    if strcmp(covname,stim.stimlegend{j})
        stimind = j;
        break;
    end
end

%is stimulus is not present set a null covops and return
if stimind==0
    %create a covops structure
    covops = struct();
    covops.isCausal = false;
    covops.isAcausal = false;
    covops.keepStim = 0;
    covops.isConst = 0;
    covops.covname = 'NOT PRESENT';
    causal_ind = 0;
    return;   
end

%decide how to pull stimulus and causalvec info
keepStim = false; %if true, then remove the stimulus and causal vec entry
causal_ind = 0; %entry in the causal vecs
if isConst==false
    if stimind<= ntdc
        %isTD = true; 
        isCausal = true;
        c_inds = find(stim.TDcausalvec(1,:)==1);
        causal_ind = c_inds(stimind);
        if stim.TDcausalvec(2,causal_ind)==1 %an acausal kernel needs the stimulus
            keepStim = true;
        end
        
    elseif stimind <= ntdc+ntdac
        %isTD = true; 
        isCausal = false;
        c_inds = find(stim.TDcausalvec(2,:)==1); %all acausal entries
        causal_ind = c_inds(stimind-ntdc);
        if stim.TDcausalvec(1,causal_ind)==1 %an causal kernel needs the stimulus
            keepStim = true;
        end
        
    end
end

%create a covops structure
covops = struct();
covops.isCausal = isCausal;
covops.isAcausal = ~isCausal;
covops.keepStim = keepStim;
covops.isConst = isConst;
covops.covname = covname;
