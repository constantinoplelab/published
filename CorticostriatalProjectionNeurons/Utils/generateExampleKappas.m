function [Kappa_MF, Kappa_MB] = generateExampleKappas(A)
%UNTITLED15 Summary of this function goes here
%   Detailed explanation goes here

rng(24) % Set random seed for reproducability

% Generate rewards
% Block structure goes mixed, low, mixed, high, mixed
r = [randsample(2:6, 40, true),...
    randsample(2:4, 40, true),... 
    randsample(2:6, 40, true),...
    randsample(4:6, 40, true),...
    randsample(2:6, 40, true)]';
r(end-20:end-15) = randsample(2:4, 6, true);

% Set up ratTrial for model to run
A_mf = A; % Preallocate A for retrospective model
A_mf.reward = 2.^r; % Add our rewards
A_mf.wait_time = ones(size(A_mf.reward));  % Set all wait times to non-nan
A_mf.prob_catch = zeros(size(r)); % Set prob catch to non-nan
A_mf.tau = 2.5; % Set tau (any non-nan value works)
A_mf.ntrials = length(r); % Only one fake session
A_mf.block = [1*ones(1,40), 3*ones(1,40), 1*ones(1,40),...
    2*ones(1,40), 1*ones(1,40)]'; % Blocks

rvec = [5 10 20 40 80]; % Reward vector
A_mb = A; % Preallocate A for inferential model
A_mb.reward = rvec(r-1)'; % Convert rewards from 1:5 to rvec
A_mb.wait_time = ones(size(A_mf.reward)); % Set all wait times to non-nan
A_mb.prob_catch = zeros(size(r)); % Set prob catch to non-nan
A_mb.tau = 2.5; % Set tau (any non-nan value works)
A_mb.ntrials = length(r); % Only one fake session

params_mf = [0.1042, 0.2527]; % Parameters for retro. simulation
params_mb = [0.6, 0.9, 0.3 0.9]; % Parameters for infer. simulation

% Generate retrospective model kappa
[~, ~, Kappa_MF] =...
    GenerateLatencyData_AlphaPHModel(params_mf, A_mf,...
    false, 'logn', nan);
Kappa_MF = Kappa_MF./max(Kappa_MF); % normalize, just for plotting purposes

% Generate inferential model kappa
[~, ~, ~, ~, ~, Kappa_MB] =...
    GenerateSynthData_Bayes(params_mb, A_mb,...
    'logn', false, nan);

end