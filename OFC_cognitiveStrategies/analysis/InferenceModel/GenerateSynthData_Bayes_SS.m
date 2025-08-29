function [WTOpt, WTMdl, BlkInf, Belief, Prior, Kappa] =...
    GenerateSynthData_Bayes_SS(params, ratTrial,...
    noise_arg, gen_arg, noise_params, rewType)
%GenerateSynthData_Bayes - generates synthetic wait times
%for optimal Bayes model. AVM, adapted by SSS
% Inputs:
%   params - Parameter vector. [kappa_mi, kappa_hi, kappa_lo, D, lambda]
%       kappa_mi, kappa_hi, kappa_lo - opportunity costs for each block.
%       D - scaling parameter 
%           ALL PARAMS ARE IN [0, 1]
%   ratTrial - Struct containing rat's data. Saved under server at
%       BehavioralModel/Data/ratTrial_ratname
%   noise_arg - What type of noise model. Always use 'logn'
%   gen_arg - Whether to generate noisy data (WTMdl) or not. Logical
%   noise_params - Variance of the noise. Always use 8
% Outputs:
%   WTOpt - Optimal wait times
%   WTMdl - Noisy wait times generate from WTOpt with constant variance
%   BlkInf - Inferred block
%   Belief - Posterior probability distribution over blocks for each trial
%   Prior - Prior probability distribution over blocks for each trial
%% Pull Parameters

H0 = 1/40;

try
    tau = ratTrial.tau;
catch
    tau = 2.5;
end

C = max(1-ratTrial.prob_catch);
min_rew = min(ratTrial.reward); 

kappa_mi = params(1)*C*log2(min_rew)/tau;
kappa_hi = params(2)*C*log2(min_rew)/tau;
kappa_lo = params(3)*C*log2(min_rew)/tau;
D = 10*params(4);

%% Predefined Probailities

p_rew = [1/5 1/5 1/5 1/5 1/5;...
    0 0 1/3 1/3 1/3;...
    1/3 1/3 1/3 0 0];

p_b = [1-H0, H0, H0;...
    H0/2, 1-H0, 0;...
    H0/2 0, 1-H0];

%% Pull Data from Rat Trial

% Only use trials where we have data from the rat (basically covers all
% non-violation trials)
isgood = ~isnan(ratTrial.wait_time);

% Pull Rewards for good trials
% Rewards are coded 1:5 - use volumes later for WT calculation
rew_vol = ratTrial.reward(isgood);

% Converts to 1:5 no matter for minvol = 4 or 5
rew_idx = convertreward(rew_vol); 

% Find when the rat starts a new sessions - restarts beliefs etc.
newday = nan(1, length(ratTrial.wait_time));
ntrials = [0; cumsum(ratTrial.ntrials)];

for ii = 1:length(ntrials)-1
    newday(ntrials(ii)+1:ntrials(ii+1)) = ii;
end

NewSessIdx = [1, find(diff(newday(isgood))) + 1];

%% Calculating Belief

[belief, prior] = deal(nan(3, length(rew_idx)));

for rr = 1:length(rew_idx)
    if ismember(rr, NewSessIdx)
        belief(:,rr) = [1/2 1/4 1/4];
    else
        prior(:,rr) = p_b*belief(:, rr-1);
        belief(:,rr) = p_rew(:, rew_idx(rr)).*prior(:,rr);
        belief(:,rr) = belief(:,rr)./sum(belief(:,rr));
    end
end

[~, blk_inf] = max(belief);

%% Calculating Wait Time

if strcmpi(rewType, 'log')
    R = log2(rew_vol);
elseif strcmpi(rewType, 'linear')
    R = rew_vol;
end

kappa_vec = [kappa_mi, kappa_hi, kappa_lo];
kappa = arrayfun(@(b) kappa_vec(b), blk_inf)';

wt_arg = (R-kappa*tau)./(kappa*tau) .* C./(1-C);
wt_opt = D*tau*log(wt_arg);

BlkInf = nan(length(ratTrial.reward), 1);
BlkInf(isgood) = blk_inf;

WTOpt = nan(length(ratTrial.reward), 1);
WTOpt(isgood) = wt_opt;

Belief = nan(3, length(ratTrial.reward));
Belief(:, isgood) = belief;

Prior = nan(3, length(ratTrial.reward));
Prior(:, isgood) = prior;

Kappa = nan(length(ratTrial.reward), 1);
Kappa(isgood) = kappa;

if gen_arg
    WTMdl = generate_wts(WTOpt, 1, noise_arg, noise_params);
else
    WTMdl = nan;
end

end