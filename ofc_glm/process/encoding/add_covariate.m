function [stim_new] = add_covariate(cnew,stim_old,covops)
%adds a covariate to an existin stimulus sctruct
%inputs: cnew: new test covariate. size ntest x nt (if timde dep), 
%                       or ntest x 1 if a constant covariate
%   stim_old: old stim structure for testing/training data
%   covops: options for chaning how this script performs
%       needs   
%               isConst: (is a stimulus or triggers a temporal kernel)
%               isCausal: (kernel goes forward in time)
%               isAcausal: (kernel goes backward in time. not mutually
%               exclusive: with isCausal. can have both as true for 2 stim)
%               covname: name of covariate/kernel
%
%
%
%outputs: stim_test,stim_train. new, updated stim structs with
%ctest_new,ctrain_n incorporated


isConst = covops.isConst; %is it a constant bias, like baseline firing rate?
%the causal and acausal parts of not exclusive
isCausal = covops.isCausal; %true, then populate correct part of causal vec
isAcausal = covops.isAcausal; %true, then populate acausal part of causal vec
covname = covops.covname; 

stim_new = stim_old; %initialize


%important sizing information
if isempty(stim_new.TDcausalvec)
    ntdc = 0;
    ntdac = 0;
else
    ntdc = sum(stim_new.TDcausalvec(1,:)); %number time-dep. causal stim
    ntdac = sum(stim_new.TDcausalvec(2,:)); %number time-dep acausal stim
end

nc = ntdc+ntdac;



%add covariates to either xmat, or TDstim. update stimlegend or xorder
if isConst %or a constant bias param
    [nx,ntrial] = size(stim_old.xmat);

    xmat = zeros(nx+1,ntrial);
    xmat(nx,:) = stim_old.xmat;
    xmat(nx+1,:) = cnew;
    stim_new.xmat = xmat;

    stim_new.Xorder{end+1} = covname;

else
    
    stimlegend_old = stim_old.stimlegend;       
        
    stim_new.TDstim{end+1} = cnew;

     if isCausal && isAcausal
         stim_new.TDcausalvec = [stim_old.TDcausalvec,[1;1]];

         %update stimlegend. one after causal. one after acausal
         stimlegend = cell(nc+2,1);
         stimlegend(1:ntdc) = stimlegend_old(1:ntdc); %add td
         stimlegend{ntdc+1} = covname; %keep original name,causal
         stimlegend(ntdc+2:ntdc+ntdac+1) = stimlegend_old(ntdc+1:ntdc+ntdac);  %add td acausal   
         stimlegend{ntdc+ntdac+2} = strcat([covname,' (ac)']); %append
         stimlegend(ntdc+ntdac+3:end) = stimlegend_old(ntdc+ntdac+1:end); %add SE %TODO: REMOVE?

     elseif isCausal && ~isAcausal
         stim_new.TDcausalvec = [stim_old.TDcausalvec,[1;0]];

         %update stimlegend. after rest of TD causal
         stimlegend = cell(nc+1,1);
         stimlegend(1:ntdc) = stimlegend_old(1:ntdc);
         stimlegend{ntdc+1} = covname;
         stimlegend(ntdc+2:end) = stimlegend_old(ntdc+1:end); %TODO: REMOVE?

     elseif ~isCausal && isAcausal
         %add 
         stim_new.TDcausalvec = [stim_old.TDcausalvec,[0;1]];

         %add stimlegend. after causal and acausal TD
         stimlegend = cell(nc+1,1);
         stimlegend(1:ntdc+ntdac) = stimlegend_old(1:ntdc+ntdac);
         stimlegend{ntdc+ntdac+1} = covname;
         stimlegend(ntdc+ntdac+2:end) = stimlegend_old(ntdc+ntdac+1:end); %TODO: REMOVE?

     end        

    stim_new.stimlegend = stimlegend;

end


    
