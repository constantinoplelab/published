function [WTOpt, WTMdl, BlkInf, RPE, Belief, Kappa, Prior] =...
    GenerateSynthData_Bayes_subopt_onlylambda_flat(params, ratTrial,...
    noise_arg, gen_arg, noise_params)
%GenerateSynthData_Bayes_subopt_onlylambda_flat - generates synthetic wait 
% times for Bayes SubOpt model with lambda over a 
% FLAT PRIOR = [1/3 1/3 1/3].
% Inputs:
%   params - Parameter vector. [kappa_mi, kappa_hi, kappa_lo, D, lambda]
%       kappa_mi, kappa_hi, kappa_lo - opportunity costs for each block.
%       D - scaling parameter 
%       lamdba - weighting between optimal prior (lambda = 1) and
%           uninformative prior (lambda = 0)
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
%   RPE - Difference between reward offer and inferred block
%   Belief - Posterior probability distribution over blocks for each trial
%   Kappa - Vector of which kappa was chosen
%   Prior - Prior probability distribution over blocks for each trial
%% Pull Parameters

H0 = 1/40;

tau = ratTrial.tau;
C = max(1-ratTrial.prob_catch);

% need min reward to normalize kappas
% But min reward is different for males (5) and females (4)
min_rew = min(ratTrial.reward); 

kappa_mi = params(1)*C*log2(min_rew)/tau;
kappa_hi = params(2)*C*log2(min_rew)/tau;
kappa_lo = params(3)*C*log2(min_rew)/tau;
D = params(4)*10;
lambda = params(5);

%% Predefined Probailities
p_rew = nan(3, 5);
p_rew(1,:) = [1/5 1/5 1/5 1/5 1/5]; % mixed block
% High and low blocks, respectively
p_rew(2,:) = [0 0 1/3 1/3 1/3];
p_rew(3,:) = [1/3 1/3 1/3 0 0];
    
p_b = [1-H0, H0, H0;...
    H0/2, 1-H0, 0;...
    H0/2 0, 1-H0];

%% Pull Data from Rat Trial

% Only use trials where we have data from the rat (basically covers all
% non-violation trials)

isgood = ~isnan(ratTrial.wait_time);
%isgood = 1:length(ratTrial.wait_time);

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

belief = nan(3, length(rew_idx));
prior = nan(3, length(rew_idx));

for rr = 1:length(rew_idx)
    if ismember(rr, NewSessIdx)
        belief(:,rr) = [1/3 1/3 1/3];
        prior(:,rr) = [1/3 1/3 1/3]; 
    else
        likey = p_rew(:, rew_idx(rr));
        prior(:,rr) = lambda*(p_b*belief(:, rr-1)) +...
            (1-lambda)*([1/3 1/3 1/3]');
        
        belief(:,rr) = likey.*prior(:,rr);
        belief(:,rr) = belief(:,rr)./sum(belief(:,rr));
    end
end

[~, blk_inf] = max(belief);

%% Calculating Wait Time

R = log2(rew_vol);

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

Kappa = nan(size(ratTrial.reward));
Kappa(isgood) = kappa;

Prior = nan(3, length(ratTrial.reward));
Prior(:, isgood) = prior;

if strcmpi(noise_arg, 'exg')
    noise_arg = 'exg_wt';
end

if gen_arg
    WTMdl = generate_wts(WTOpt, 3, noise_arg, noise_params);
else
    WTMdl = nan;
end

%% Calculating RPE

RPE = nan(length(ratTrial.reward), 1);
RPE(isgood) = rew_vol - kappa;

end
