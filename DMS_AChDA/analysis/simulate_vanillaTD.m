function [LatOpt, RPE, Kappa, Rewards, Delays] =...
    simulate_vanillaTD(alpha, ratTrial, rpetype, scalingfactor)

if nargin<4
    scalingfactor = 0.48;
    fprintf('Using default discount factor %.2f\n', scalingfactor)
end
%% Pull Data from Rat Trial

% Only use trials where we have data from the rat (basically covers all
% non-violation trials)
isgood = ~isnan(ratTrial.ITI);

Rewards = nan(size(ratTrial.reward));
Rewards(isgood) = ratTrial.reward(isgood);

Delays = nan(size(ratTrial.reward_delay));
Delays(isgood) = ratTrial.reward_delay(isgood);

% Pull Rewards/Delays for good trials
% Rewards are coded 1:5 - use volumes later for ITI calculation
if strcmpi(rpetype, 'reward')
    newrew = ratTrial.reward(isgood);
elseif strcmpi(rpetype, 'delay')
    newrew = ratTrial.reward_delay;
    newrew(ratTrial.catch==1) = 0;
    newrew = newrew(isgood);
end

% Find when the rat starts a new sessions - restarts beliefs etc.
newday = nan(1, length(ratTrial.ITI));
ntrials = [0; cumsum(ratTrial.ntrials)];

for ii = 1:length(ntrials)-1
    newday(ntrials(ii)+1:ntrials(ii+1)) = ii;
end

NewSessIdx = [1, find(diff(newday(isgood))) + 1];

%% Calculating Belief

% Initalize vectors
kappa = nan(length(newrew), 1);
rpe = nan(length(newrew), 1);

% Initial conditions
kappa_initial = 1;

%%run a first pass to determine the initial values of kappa
for rr = 1:length(newrew)-1
    % First trial - restart
    if ismember(rr, NewSessIdx)
        kappa(rr) = kappa_initial;
    end
    rpe(rr) = newrew(rr) - kappa(rr);
    kappa(rr+1) = kappa(rr) + alpha * rpe(rr);
end

%%now, use the distribution of kappa values to choose a reasonable initial
%%condition. Then rerun to get the new kappas/rpeRs
kappa_initial = median(kappa); 
kappa = nan(length(newrew), 1);
rpe = nan(length(newrew), 1);
for rr = 1:length(newrew)-1
    % First trial - restart
    if ismember(rr, NewSessIdx)
        kappa(rr) = kappa_initial;
    end

    rpe(rr) = newrew(rr) - kappa(rr);
    kappa(rr+1) = kappa(rr) + alpha * rpe(rr);
end

%% Calculating ITI

lat_opt = prctile(ratTrial.ITI, 90) - scalingfactor*kappa;  

LatOpt = nan(length(ratTrial.reward), 1);
LatOpt(isgood) = lat_opt;

RPE = nan(length(ratTrial.reward), 1);
RPE(isgood) = rpe;

Kappa = nan(size(ratTrial.reward));
Kappa(isgood) = kappa;


end