function newtrials = balanceTrials(mask1,mask2,trials)
    
    %balances trials w.r.t. a binary mask.
    %i.e., if left and right trials need to be balanced
    %mask1 gives true for all left trials, and false otherwise
    %this will give new data set with equal amounts of each, via
    %subsampling
    
    %mask1, mask2, and trials should be the same length
    
    inds1 = find(mask1==1);
    inds2 = find(mask2==1);
    n1 =numel(inds1);
    n2 = numel(inds2); 
    
    if n1 >n2
        inds1 = inds1(1:n2);
        newtrials = false(size(trials));
        newtrials(inds1) = true;
        newtrials(inds2) = true;
    elseif n2 > n1
        inds2 = inds2(1:n1);
        newtrials = false(size(trials));
        newtrials(inds1) = true;
        newtrials(inds2) = true;
    else
        newtrials = trials;
    end
       
end