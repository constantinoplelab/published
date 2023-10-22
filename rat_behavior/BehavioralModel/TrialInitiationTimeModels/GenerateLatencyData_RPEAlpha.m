function [LatOpt, LatMdl, RPE, Kappa] =...
    GenerateLatencyData_RPEAlpha(params, ratTrial,...
    gen_arg, noise_arg, noise_params)
%GenerateLatencyData_VanillaAlpha
%   params - model parameters
%       params(1) = learning rate, in [0, 1]
%       params(2) = scale parameter, [0, 1]
%   ratTrial - trial data. Usually loaded from ProcessedRatData/A_Structs
%   gen_arg - if you want to generate noisy data. If false, LatMdl = nan.
%   noise_arg - noise model. Use 'logn'.
%   noise_params - the value of the noise variance. Use 4.
%% Pull Parameters

alpha0 = params(1);
D = params(2)*15;

%% Pull Data from Rat Trial

% Only use trials where we have data from the rat (basically covers all
% non-violation trials)
isgood = ~isnan(ratTrial.ITI);

% Pull Rewards for good trials
% Rewards are coded 1:5 - use volumes later for WT calculation
rew = ratTrial.reward(isgood);
rew = log2(rew);

% Find when the rat starts a new sessions - restarts beliefs etc.
newday = nan(1, length(ratTrial.ITI));
ntrials = [0; cumsum(ratTrial.ntrials)];

for ii = 1:length(ntrials)-1
    newday(ntrials(ii)+1:ntrials(ii+1)) = ii;
end

NewSessIdx = [1, find(diff(newday(isgood))) + 1];

%% Calculating Belief

% Initalize vectors
kappa = nan(length(rew), 1);
rpe = nan(length(rew), 1);
alpha = nan(length(rew), 1);

% Initial conditions
kappa(1) = 0.5;

for rr = 1:length(rew)-1
    % First trial - restart
    if ismember(rr, NewSessIdx)
        kappa(rr) = kappa(1);
    end

    rpe(rr) = rew(rr) - kappa(rr);
    alpha(rr) = abs(rpe(rr))*alpha0;
    
    if alpha(rr) > 1
        alpha(rr) = 1;
    end
    
    kappa(rr+1) = kappa(rr) +  alpha(rr) * rpe(rr);
end

%% Calculating ITI
lat_opt = D./kappa;

LatOpt = nan(length(ratTrial.reward), 1);
LatOpt(isgood) = lat_opt;

RPE = nan(length(ratTrial.reward), 1);
RPE(isgood) = rpe;

Kappa = nan(size(ratTrial.reward));
Kappa(isgood) = kappa;

if strcmpi(noise_arg, 'exg')
    noise_arg = 'exg_iti';
end

if gen_arg
    LatMdl = generate_wts(LatOpt, 1, noise_arg, noise_params);
else
    LatMdl = nan;
end

end
