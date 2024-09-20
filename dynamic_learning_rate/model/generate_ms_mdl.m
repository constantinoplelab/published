function [itiMdl, itiOpt, V, rpe, alpha, G] = generate_ms_mdl(...
    R, alpha0, D)
% generate_ms_mdl - generates data for Mackintosh Surprise model
% INPUTS:
%   R - reward vector. Must be 1:5
%   alpha0 - Base learning rate
%   D - Scale parameter
% OUTPUTS:
%   itiMdl - trial intiation time with added noise
%   itiOpt - optimal trial intiation time
%   V - value estimate for each trial
%   rpe - RPE for each trial
%   alpha - Learning rate for each trial
%   G - Gain on learing rate for each trial

% Preallocate data strucutres
[V, rpe, alpha] = deal(nan(size(R)));

% Intial conditions
V(1) = 2.5;

alphas = alpha0*(1:5); % alpha proprotional to reward offer

% Loop over trials
for tt = 1:length(R)
    alpha(tt) = alphas(R(tt));

    rpe(tt) = R(tt) - V(tt);
    V(tt+1) = V(tt) + alpha(tt)*rpe(tt);
end
G = alpha./alpha0;

itiOpt = D./V; % Optimal ITI
itiMdl = generate_wts(itiOpt, 1, 'logn', 1.5); % ITI with noise

end

