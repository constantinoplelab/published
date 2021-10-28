function [trial_idx,idx_flat] = parse_pseudotrials(ops,trial_id)
%make pseudotrials of continuous-time trial sets
%based on trial_idx binary mask, will look at start and
%stop indices in ops and 
%
%trial_id: trial numbers to be used/extracted from data
%
%output
%trial_idx: [ns,2] array of ns "pseudotrials" of continous data. 
%   trial_idx[:,1] is start index of that trial, trial_idx[:,2] is end
%   index
%
%idx_flat: flattened version of trial_idx, giving index of every data point


%find the start,end times from ops, break entire sequence into "trials"
trial_id = ops.Dat.trial_idx(trial_id); %convert from logical to integer
start_idx = ops.Dat.start_idx(trial_id);
stop_idx= ops.Dat.end_idx(trial_id);

%reshape to ensure column vectors
trial_id = reshape(trial_id,numel(trial_id),[]);
start_idx = reshape(start_idx,numel(start_idx),[]);
stop_idx = reshape(stop_idx,numel(stop_idx),[]);

%find contiguous trials and reform into a pseudo-trial structure
%with continuous trials acting as one trial, and breaks with time
%breaks
trialbreaks = find(diff(trial_id) > 1);
      
if ~isempty(trialbreaks)
    %concat start_idx of continous trial with stop_idx of continous trial
    trial_idx = [start_idx(1),stop_idx(trialbreaks(1));
                start_idx(trialbreaks(1:end-1)+1),stop_idx(trialbreaks(2:end));
                start_idx(trialbreaks(end)+1),stop_idx(end)];
            
else %case of no trial breaks.
    tstart = start_idx(1);   
    %tend = (stop_idx(end)-start_idx(1))+1; %+1?   
    tend = stop_idx(end) ;
    trial_idx = [tstart,tend];

end

%make flattened version to extract only relevant data from stimuli, spike
%times, etc.
ns = size(trial_idx,1);
ntns = sum(diff(trial_idx,1,2)+1);
idx_flat = zeros(ntns,1);
ind = 0;
for j = 1:ns
    nt = numel(trial_idx(j,1):trial_idx(j,2));
    idx_flat(ind + 1:ind + nt) = trial_idx(j,1):trial_idx(j,2);
    
    ind = ind + nt;
end



        