function [normval, rbar] = divisiveNorm_fun_2(A, params, rewType)
% Divisive normalization model. Based on model in Khaw et al. 2017 https://doi.org/10.1073/pnas.1715293114

%INPUTS:
%   A = .mat file of rat behavior data
%   params = 3 parameter values (K, alpha, N)
%   rewType = string indicating whether to use log2 rewards or linear
%       rewards

%% Parameters
k = params(1);
alpha = params(2); 
tback = params(3); %number of trials to sum over (N)

%% Pull rat data

% Pull Rewards
rew_vol = A.reward;

% Determine if rat has 4 or 5ul rewards, will be used to set intial sum
if rem(rew_vol(end), 5) == 0
    trialAvg = log2(30.6); %empirically determined average offer per trial for males
else 
    trialAvg = log2(24.7); %empirically determined average offer per trial for females
end

% Find when the rat starts a new sessions - restarts trial integration
% window
newday = nan(1, length(A.wait_time));
ntrials = [0; cumsum(A.ntrials)];

for ii = 1:length(ntrials)-1
    newday(ntrials(ii)+1:ntrials(ii+1)) = ii;
end

NewSessIdx = [1, find(diff(newday)) + 1];

%% calculate normalized value 
vt = nan(length(rew_vol), 1);
rbar = vt;
normval = vt; %wait time is proportional to normalized value

if strcmpi(rewType, 'log')
    vt = log2(rew_vol); %consistent with Mah 2023, assumes rats use commpressed utility functions
elseif strcmpi(rewType, 'linear')
    vt = rew_vol;
end

for j = 1:length(rew_vol)
     if ismember(j, NewSessIdx) 
        t = 1; %reset at the start of each new session
        idx = j;
     else 
        t = t+1;
     end

    if t>tback
        rbar(j) = sum(vt(j-tback:j), 'omitnan');
    else %reset to average and then start summing new offers
        rbar(j) = (tback-t)*trialAvg + sum(vt(idx:j), 'omitnan'); %includes trial average in the window until t = tback
    end
    normval(j) = k*(vt(j)./(1+alpha*rbar(j))); %divisive normalization
end
