function [wtNorm_mean, wtRaw_mean] = prev10Rews(A, ntrials, inferredMixArg)
% Compute wait times conditioned on the sum over the previous (ntrials)
% trials

%INPUTS: 
%   A = .mat file of individual rat behavior data
%   ntrials = number of trials to sum over
%   inferredMixArg = whether to also filter by trials inferred to be in a mixed
%       block by a bayesian inference agent

%OUTPUTS:
%   wtNorm_mean = normalized wait times for 20ul trialsconditioned on the
%       sum of the previous trials

[~, rew] = convertreward(A.reward);
wt = A.wait_time;

trials = 1:length(A.block);

% Even the bayes inference model occassionally makes mistaken inferences 
% in mixed blocks after long strings of only low or high volume trials.
% Do not include trials that are inferred to be in a high or low block. 
% Note: this likely does not perfectly capture mistaken inferences by each
% rat
if ~isempty(inferredMixArg)
    BlkInf = inferredMixArg;
    find_20ul_catch = find(rew==20 & A.block==1 &...
        A.optout==1 & A.catch==1 & BlkInf==1); %20ul catch trials that are correctly inferred to be mixed block trials based on inference model
else
    find_20ul_catch = find(rew==20 & A.block==1 &...
        A.optout==1 & A.catch==1);
end

find_20ul_catch(find_20ul_catch < ntrials) = []; %don't include catch trials before trial window to look at

% get indices of the last 10 trials before the 20ul catch trial
prevHits = arrayfun(@(x) trials(find(trials < find_20ul_catch(x), ntrials, 'last')), ...
    1:length(find_20ul_catch), 'UniformOutput', false); 

% don't include 20ul trials that don't have (ntrials) previous trials
bad = arrayfun(@(x) length(prevHits{x}) < ntrials, 1:length(prevHits));

prevHits(bad) = [];
find_20ul_catch(bad) = [];

% compute sum of previous trials for each 20ul trial
prevRew = arrayfun(@(x) sum(rew(prevHits{x})), 1:length(prevHits));

%find bottom and top 50% of sums
p_low = prctile(prevRew, 50);
p_high = prctile(prevRew, 50);

less = prevRew < p_low;
greater = prevRew > p_high;

%Wait times for bottom and top 50% of sums
wtNorm = (wt - mean(wt(find_20ul_catch), 'omitnan'))./...
    std(wt(find_20ul_catch), 'omitnan');

wtNorm_mean = [mean(wtNorm(find_20ul_catch(less)), 'omitnan'),...
    mean(wtNorm(find_20ul_catch(greater)), 'omitnan')];

wtRaw_mean = [mean(wt(find_20ul_catch(less)), 'omitnan'),...
    mean(wt(find_20ul_catch(greater)), 'omitnan')];
