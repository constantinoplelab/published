function [Dat,stim,ops] = parse_stimuli_wrapper(ops)
%decides which parse_stimuli function to use. alters stim for CPD calc.
%
% inputs: ops.stimwrapper decides which function to use
%
%

switch ops.stimwrapper
case 'ofc'
    %[stim_train,stim_test] = parse_stimuli_ofc(ops);
    if ops.usePseudotrial
        [Dat,stim] = loadGlmDat_ofc(ops,ops.dataname);
    else
        [Dat,stim] = loadGlmDat_ofc_trialbased(ops,ops.dataname);
    end

case 'gt' 
    stim = parse_stimuli_gt(ops);
    Dat = ops.Dat; 
end

%decide on any post-hoc changes, such as omittied a kernel, or omitting
%behavior
%decide if any post-hoc changes to stimulus should be made

%omit a kernel
if isfield(ops,'omitKernels') && numel(ops.omitKernels) > 0
    dispL('omitting kernels from regular stimtype',3,ops.dispL);
    n_omit = numel(ops.omitKernels);
    for o_ind = 1:n_omit
        stim = remove_covariate(stim,ops.omitKernels{o_ind});
    end
end

%remove behavior
if isfield(ops,'omitBehavior') && numel(ops.omitBehavior) > 0
    %dispL('removing behavior from kernels by merging two kernels',3,ops.dispL);
    n_merge = size(ops.omitBehavior,1);
    for o_merge = 1:n_merge
        
        %find number of covariates to merge, and which ones are new
        nkern_check = size(ops.omitBehavior(o_merge,:),2);
        nold = 0;
        causal_ind_vec = zeros(nkern_check,1);
        for j = 1:nkern_check
            kern_j = ops.omitBehavior{o_merge,j};
            [~,caus_ind,covops] = stimulusInfo(stim,kern_j);
            if caus_ind > 0 %found in stim
                nold = nold + 1;
                causal_ind_vec(j) = caus_ind;
            end
        end
        
        %handle adding and removing causal and causal bits in one swoop:
        %if either (a)causal kernels are part of this list and have
        %a counterpart, then they are automatically removed/added.
        kern_j = ops.omitBehavior{o_merge,1};
        [~,~,covops] = stimulusInfo(stim,kern_j);
        if covops.keepStim
            covops.isCausal = 1;
            covops.isAcausal = 1;
        end
        
        %which stimulus indices should be merged. unique() takse care of
        %acausal redundancy
        caus_inds = unique(causal_ind_vec(causal_ind_vec>0));
        mergestim = zeros(size(stim.TDstim{caus_inds(1)}));
        for j = 1:numel(caus_inds)
            mergestim = mergestim + stim.TDstim{caus_inds(j)};
        end
        
        %change to new name, provided by last index in omit behavior
        covops.covname = ops.omitBehavior{o_merge,nold+1};
        %add the new covariate
        stim = add_covariate(mergestim,stim,covops);
        
        
        %remove the old ones
        for j = 1:nold
            kern_j = ops.omitBehavior{o_merge,j};
            stim = remove_covariate(stim,kern_j);
        end
        
        %{
        kern1 = ops.omitBehavior{o_merge,1};
        kern2 = ops.omitBehavior{o_merge,2};
        [~,causal_ind1,covops1] = stimulusInfo(stim,kern1);
        [~,causal_ind2,~] = stimulusInfo(stim,kern2);

        %merged stim
        mergestim = stim.TDstim{causal_ind1} + stim.TDstim{causal_ind2};

        %change to new name, provided by 3rd index in omit behavior
        covops1.covname = ops.omitBehavior{o_merge,3};
        %add the new covariate
        stim = add_covariate(mergestim,stim,covops1);
        %remove the old ones
        stim = remove_covariate(stim,kern1);
        stim = remove_covariate(stim,kern2);
        %}

    end
end