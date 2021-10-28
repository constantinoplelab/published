function [ordering,groupbreaks] = sortAmplitudeByCluster(X,grouping)
%does time 2 peak by subgroups of data
%
%Inputs:
%X: (nn x nt) data of interest, nn neurons, nt time points
%ops: structure with following
%   useZScore: do normal sorting on X, or z score it first
%   grouping: (nn x 1) vector of integers, assignign each to a class
%
%Outputs:
%ordering: (nn x 1) vector of idx of each neurons, sorted by time 2 peak,
%but within in each cluster

%find cluster boundaries in grouping
nn = numel(grouping);
[~,I_order] = sort(grouping);
groupbreaks = find(diff(sort(grouping))); %when does a new group start
ordering = zeros(nn,1);


%loop over groups and find time to peak
changevec_temp = [1,groupbreaks,nn+1];
inds = 0;
for m = 1:numel(groupbreaks)+1
    %neuron_idx in that group
    neur_idx = I_order(changevec_temp(m):changevec_temp(m+1)-1); 
    n_neur = numel(neur_idx);

    %ordering_m = time2peak(X(neur_idx,:),ops);
    [~,ordering_m] = sort(X(neur_idx,:),'descend');
    ordering(inds+1:inds+n_neur) = neur_idx(ordering_m);

    inds = inds+n_neur;
end

