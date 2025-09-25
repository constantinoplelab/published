function [modelITI, modelValue, ratTrial, modelRPE] = ...
    vanillaTD_perRat(A, modeltype, alpha, scalingfactor)
% Simulate trial initiation time by block using a vanilla-alpha TD
% algorithm.

if nargin<4
    scalingfactor = 0.48;
    fprintf('Using default scaling factor %.2f\n', scalingfactor)
    if nargin<3
        alpha = 0.5;
        fprintf('Using default learning rate %.3f\n', alpha)
    end
end

ratTrial = A;

% convert reward from uL to nominal scale
[ratTrial.reward,~] = convertreward(ratTrial.reward);

% Process trial initiation time data 
% set trial initiation on first trial of each session to NaN
ctr = 1;
for jk = 1:length(ratTrial.ntrials)
    ratTrial.ITI(ctr) = nan;
    ctr = ctr+ratTrial.ntrials(jk);
end
% remove outliers
disp('Removing bottom,top 1% in ITI..')
thresh = [prctile(ratTrial.ITI,1), prctile(ratTrial.ITI,99)];
ratTrial.ITI(ratTrial.ITI<min(thresh) | ratTrial.ITI>max(thresh)) = nan;

% Simulate trial initiation time with a vanilla TD model
[modelITI, modelRPE, modelValue, ~, ~] =...
    simulate_vanillaTD(alpha, ratTrial, modeltype, scalingfactor);


end

