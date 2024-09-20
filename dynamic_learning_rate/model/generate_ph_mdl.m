function [itiMdl, itiOpt, V, rpe, alpha, G] = generate_ph_mdl(...
    R, alpha0, D)
% generate_ph_mdl - generates data for Pearce-Hall model
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

% Preallocate data structures
[V, rpe, alpha, G] = deal(nan(size(R)));

% Intial condition
V(1) = 2.5;

% Loop over trials
for tt = 1:length(R)
    if tt == 1
        G(tt) = 1;
    else
        G(tt) = abs(rpe(tt-1)); % Gain is unsigned previous RPE
    end
    alpha(tt) = min(alpha0*G(tt), 1);

    rpe(tt) = R(tt) - V(tt);
    V(tt+1) = V(tt) + alpha(tt)*rpe(tt);
end

itiOpt = D./V; % Optimal ITI
itiMdl = generate_wts(itiOpt, 1, 'logn', 1.5); % ITI with noise

end