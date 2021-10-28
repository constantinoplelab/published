function kfold_inds = cvpartition_wrapper(ops)
%decides how to generate labels for cross-validation
%
%



if ops.kfold > 1
    switch ops.cvtype
        case 'ofc'
            %use CV partition, but fit it into your current code
            %decide on stratification groups. currently balance among trial behaviors
            ofcdat = load(ops.dataname);
            group = cvalgroup(ofcdat,'behavior');
            %group = cvalgroup(ofcdatat,'prob');
            
            if ops.removeData %removes opt-out trials
                group = group(ops.Dat.trial_idx);
            end

            %create object with a group
            rng('default'); %otherwise cvpartition isn't reproducible. has to go RIGHT before cvpartition call
            kfold_obj = cvpartition(group,'KFold',ops.kfold);

        case 'default' %default case of splitting by trial number
            ns = numel(ops.Dat.trial_idx);
            rng('default'); %otherwise cvpartition isn't reproducible. has to go RIGHT before cvpartition call
            kfold_obj = cvpartition(ns,'KFold',ops.kfold);
            
        case 'continuous' %cross-validated sets that are trial-continuous
            ns = numel(ops.Dat.trial_idx);
            kfold_obj = struct();
            indsvec = zeros(ns,1);
            indsvec(1:ns < ns/ops.kfold) = 1;
            for j = 2:ops.kfold-1
                indsvec(1:ns <= (j)*ns/ops.kfold & 1:ns >= (j-1)*ns/ops.kfold) = j;
            end
            indsvec(1:ns > (ops.kfold-1)*ns/ops.kfold) = ops.kfold;
            
            testinds = @(k) indsvec == k;
            
            kfold_obj.test = testinds;
                        
            
        case 'stimbalance' %test case to balance stimuli instances
            nkern = numel(ops.Dat.covariates);
            ns = numel(ops.Dat.trial_idx);
            group_base = false(ns,nkern);
            for j = 1:ns
                for k = 1:nkern
                    if find(ops.Dat.covariates{k}(:,j)) > 0
                        group_base(j,k) = true;
                    end
                end
            end
            group = bi2de(group_base);
        
            rng('default'); %otherwise cvpartition isn't reproducible. has to go RIGHT before cvpartition call
            kfold_obj = cvpartition(group,'KFold',ops.kfold);
            

    end
    
    ns = numel(ops.Dat.trial_idx);
    kfold_inds = zeros(ns,1);
    for j = 1:ops.kfold
        kfold_inds(kfold_obj.test(j))=j;
    end
    
else  %no paritioning. probably not advisable, but good test condition
    ns = numel(ops.Dat.trial_idx);
    kfold_inds = ones(ns,1);
end