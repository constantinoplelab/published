function RunLengthPosterior = calculate_changepoint_model(R, Ntruc)
% RunLengthPosterior - calculate Bayesian Online Changepoint Detection
% model (Adams and McKay, 2007)
% INPUTS:
%   R - Reward vector
%   Ntruc - number of trials to trucate runlenght posterior (~70 is
%   efficient)
% OUTPUTS:
%   RunLengthPosterior - Run-length posterior on each trial (each trial is
%   a column)

H0 = 1/40; % Hazard rate

% Preallocate run-length posterior
run_length_posterior = cell(size(R));
run_length_posterior{1} = 1;

% Loop over trials
for t = 2:length(R)
    fprintf('%d out of %d\n', t, length(R)) % Display progress

    %calculates pi (step 3)
    Pi = arrayfun(@(l) calc_pi(l, t, R), 1:min(Ntruc, t-1));

    % calculate growth probabilities
    GrowthP = run_length_posterior{t-1} .* Pi * (1-H0);

    % calculate change point
    ChangePointP = sum(run_length_posterior{t-1} .* Pi .* H0);
    % Store PXandR
    PXandR = [ChangePointP, GrowthP(1:min(t-1, Ntruc-1))];

    % Normalize posterior
    run_length_posterior{t} = PXandR./sum(PXandR);
end

% Puts everything in rectangular matrix
RunLengthPosterior = nan(Ntruc, length(R));

for t = 1:length(R)
    RunLengthPosterior(1:min(t, Ntruc), t) = run_length_posterior{t};
end

end


function [Pi, PXgivenB] = calc_pi(l, t, R)

lik = [1/5 1/5 1/5 1/5 1/5;...
    0 0 1/3 1/3 1/3;...
    1/3 1/3 1/3 0 0];
pb = [1/3 1/3 1/3];

PXgivenB =...
    arrayfun(@(b) prod(arrayfun(@(tt) lik(b, R(tt)),...
    (t-l):t-1)), 1:3);
PXgivenB = PXgivenB./sum(PXgivenB);

Pi = (PXgivenB.*pb)*lik(:, R(t));

end
