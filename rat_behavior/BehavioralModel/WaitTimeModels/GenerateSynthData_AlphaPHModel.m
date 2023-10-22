function [WTOpt, WTMdl, Kappa, RPE, Belief, Alphas] =...
    GenerateSynthData_AlphaPHModel(params, ratTrial,...
    noise_arg, gen_arg, noise_params)
%GenerateSynthData_AlphaPHModel - generates synthetic wait times for 
% dynamic learning rate model
% Inputs:
%   params - Parameter vector. [alpha0, D]
%       alpha0 - baseline learning rate.
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
%   Kappa - Recency-weighted running average of preivous rewards
%   RPE - Difference between reward offer and kappa
%   Belief - Posterior probability distribution over blocks for each trial.
%       Calculated using Bayes optimal model for this case.
%   Alphas - Dynamic alpha
%% Pull Parameters

H0 = 1/40;

alpha0 = params(1); % model-free learning rate
D = 15*params(2); % scale parameter

%% Predefined Probailities
% Predefined probabilities
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

% Pull Rewards for good trials
% Rewards are coded 1:5 - use volumes later for WT calculation
rew_vol = ratTrial.reward(isgood);
rew_idx = convertreward(rew_vol); % Converts to 1:5

rew_vol = log2(rew_vol);

% Find when the rat starts a new sessions - resets beliefs etc.
newday = nan(1, length(ratTrial.wait_time));
ntrials = [0; cumsum(ratTrial.ntrials)];

for ii = 1:length(ntrials)-1
    newday(ntrials(ii)+1:ntrials(ii+1)) = ii;
end

NewSessIdx = [1, find(diff(newday(isgood))) + 1];

%% Calculating Belief

belief = nan(3, length(rew_vol));
belief(:,1) = [1/3 1/3 1/3];

rpe = nan(size(rew_vol));
alphas = nan(size(rew_vol));

kappa = nan(size(rew_vol));
kappa(1) = mean(log([5 10 20 40 80]));

for rr = 1:length(rew_vol)-1
    if ismember(rr, NewSessIdx)
        belief(:,rr) = [1/2 1/4 1/4];
        kappa(rr) = kappa(1);
        alphas(rr) = alpha0;
    else
        likey = p_rew(:, rew_idx(rr));
        prior = p_b*belief(:, rr-1);
                
        belief(:,rr) = likey.*prior;
        belief(:,rr) = belief(:,rr)./sum(belief(:,rr)); 
        alphas(rr) =...
            min(alpha0 * 1/(1-abs(belief(1,rr) - belief(1,rr-1))),...
            1);
    end
    
    rpe(rr) = rew_vol(rr) - kappa(rr);
    kappa(rr+1) = kappa(rr) + alphas(rr)*rpe(rr);
end

%% Calculating Wait Time

tau = ratTrial.tau;
C = max(1-ratTrial.prob_catch);

kappa = (kappa/max(kappa)).*(min(rew_vol)*C/tau);
kappa(kappa<0) = eps;

wt_arg = (rew_vol-kappa*tau)./(kappa*tau) .* C./(1-C);
wt_opt = D*tau*log(wt_arg);
wt_opt(wt_opt<0) = eps;

WTOpt = nan(length(ratTrial.reward), 1);
WTOpt(isgood) = wt_opt;

Belief = nan(3, length(ratTrial.reward));
Belief(:, isgood) = belief;

Kappa = nan(length(ratTrial.reward), 1);
Kappa(isgood) = kappa;

RPE = nan(length(ratTrial.reward), 1);
RPE(isgood) = rpe;

Alphas = nan(length(ratTrial.reward), 1);
Alphas(isgood) = alphas;

if gen_arg
    WTMdl = generate_wts(WTOpt, 1, noise_arg, noise_params);
else
    WTMdl = nan;
end

end