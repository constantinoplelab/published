function [ITIOpt, ITIMdl, Kappa, RPE, Belief, Alphas] =...
    GenerateLatencyData_AlphaPHModel(params, ratTrial,...
    gen_arg, noise_arg, noise_params)
%GenerateSynthData_BeliefState Generates wait times for belief state model
%% Pull Parameters

H0 = 1/40;

alpha0 = params(1); % model-free learning rate
D = params(2)*25; % scale parameter

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
        
%         surprise = calculateSurprise(prior, belief(:,rr));
%         alphas(rr) = alpha0 * surprise;
        alphas(rr) = alpha0 * 1/(1-abs(belief(1,rr) - belief(1,rr-1)));
        
        % Prevent alpha blow up
        if alphas(rr) > 1
            alphas(rr) = 1;
        end 
    end
    
    rpe(rr) = rew_vol(rr) - kappa(rr);
    kappa(rr+1) = kappa(rr) + alphas(rr)*rpe(rr);
end

%% Calculating Wait Time

kappa(kappa<0) = eps;

iti = D./kappa;

ITIOpt = nan(size(ratTrial.reward));
ITIOpt(isgood) = iti;

Belief = nan(3, length(ratTrial.reward));
Belief(:, isgood) = belief;

Kappa = nan(length(ratTrial.reward), 1);
Kappa(isgood) = kappa;

RPE = nan(length(ratTrial.reward), 1);
RPE(isgood) = rpe;

Alphas = nan(length(ratTrial.reward), 1);
Alphas(isgood) = alphas;

if gen_arg
    ITIMdl = generate_wts(ITIOpt, 1, noise_arg, noise_params);
else
    ITIMdl = nan;
end
end