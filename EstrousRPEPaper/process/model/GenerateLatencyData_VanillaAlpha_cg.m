function [LatOpt, RPE, Kappa, Rewards] =...
    GenerateLatencyData_VanillaAlpha_cg(params, ratTrial, opto_arg,...
    gain_arg, fig)

%% Pull Parameters
alpha = params(1);
if gain_arg
    rpe_gain = params(2);
end
%% Pull Data from Rat Trial

% Only use trials where we have data from the rat (basically covers all
% non-violation trials)
isgood = ~isnan(ratTrial.ITI);

% Pull Rewards for good trials
% Rewards are coded 1:5 - use volumes later for ITI calculation
newrew = ratTrial.reward(isgood);
% violations = ratTrial.vios==1;
% newrew(violations) = 0; %set rewards on violation trials to 0

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
%%condition. Then rerun to get the new kappas/rpes
kappa_initial = median(kappa); 
maxRPE = max(rpe);
largeRPEs_threshold = prctile(rpe, 90);
kappa = nan(length(newrew), 1);
rpe = nan(length(newrew), 1);
for rr = 1:length(newrew)-1
    % First trial - restart
    if ismember(rr, NewSessIdx)
        kappa(rr) = kappa_initial;
    end

    if opto_arg == 1 && rand<=0.3
        rpe(rr) = newrew(rr) - kappa(rr) + maxRPE;
    elseif gain_arg==1
        if (newrew(rr) - kappa(rr)) > largeRPEs_threshold
            rpe(rr) = (newrew(rr) - kappa(rr)) * rpe_gain;
        else
            rpe(rr) = newrew(rr) - kappa(rr);
        end
    else
        rpe(rr) = newrew(rr) - kappa(rr);
    end
    kappa(rr+1) = kappa(rr) + alpha * rpe(rr);
    
end

%% Calculating ITI

if fig==3
    lat_opt = prctile(ratTrial.ITI, 90)-0.1*kappa;  
elseif fig==4
    lat_opt = prctile(ratTrial.ITI, 84)-0.08*kappa; 
end
% lat_opt(lat_opt<2) = 2; %floor effect

LatOpt = nan(length(ratTrial.reward), 1);
LatOpt(isgood) = lat_opt;

RPE = nan(length(ratTrial.reward), 1);
RPE(isgood) = rpe;

Kappa = nan(size(ratTrial.reward));
Kappa(isgood) = kappa;

Rewards = nan(size(ratTrial.reward));
Rewards(isgood) = newrew;
end