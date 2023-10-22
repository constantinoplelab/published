function [LatOpt, LatMdl, BlkInf, RPE, Belief, Kappa] =...
    GenerateLatencyData_Bayes_subopt_onlylambda(params, ratTrial,...
    gen_arg, noise_arg, noise_params)
%%Similar to GenerateSynthData_Bayes, but with two new parameters: lambda
%%(prior weighting) and lapse rate. cmc 8/19/21.
%% Pull Parameters

H0 = 1/40;

tau = ratTrial.tau;
C = max(1-ratTrial.prob_catch);

kappa_mi = params(1)*C*log2(5)/tau;
kappa_hi = params(2)*C*log2(5)/tau;
kappa_lo = params(3)*C*log2(5)/tau;
D = params(4)*10;
lambda = params(5);

%% Predefined Probailities
% Predefined probabilities
% p_rew = [1/5 1/5 1/5 1/5 1/5;...
%     0 0 1/3 1/3 1/3;...
%     1/3 1/3 1/3 0 0];

p_rew = nan(3, 5);
p_rew(1,:) = [1/5 1/5 1/5 1/5 1/5]; % mixed block
p_rew(2,:) = [0 0 1/3 1/3 1/3]; % high block
p_rew(3,:) = [1/3 1/3 1/3 0 0]; % low block
    
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

for rr = 1:length(rew_idx)
    if ismember(rr, NewSessIdx)
        belief(:,rr) = [1/3 1/3 1/3];
    else
        likey = p_rew(:, rew_idx(rr));
        prior = lambda*(p_b*belief(:, rr-1)) + (1-lambda)*([1/3 1/3 1/3]');
        
        % unsure of what lambda should be doing... if it's weighting null 
        % prior and true prior, then I would expect something like like 72, 
        % as opposed to like 70, which is just changing the belief on the
        % previous trial.
        
        belief(:,rr) = likey.*prior;
        belief(:,rr) = belief(:,rr)./sum(belief(:,rr));
    end
end

[~, blk_inf] = max(belief);

%% Calculating Wait Time

kappa_vec = [kappa_mi, kappa_hi, kappa_lo];
kappa = arrayfun(@(b) kappa_vec(b), blk_inf)';

lat_opt = D./kappa;

% wt_arg = (R-kappa*tau)./(kappa*tau) .* C./(1-C);
% wt_opt = D*tau*log(wt_arg);

BlkInf = nan(length(ratTrial.reward), 1);
BlkInf(isgood) = blk_inf;

LatOpt = nan(length(ratTrial.reward), 1);
LatOpt(isgood) = lat_opt;

Belief = nan(3, length(ratTrial.reward));
Belief(:, isgood) = belief;

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

%% Calculating RPE

RPE = nan(length(ratTrial.reward), 1);
RPE(isgood) = rew_vol - kappa;

end
