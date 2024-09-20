
function [nLL, nLL_Trial] =...
    lognNLL_latency_secondhalf(params, lat2fit, ratTrial,...
    mdl, likparams)
%lognNLL_latency_secondhalf - calculate nLL of model parameters given data
% INPUTS:
%   params - model parameters
%   lat2fit - data to compare model outputs to
%   ratTrial - rat beahvioral data
%   mdl - model name. Use 'VanillaAlpha'
%   likparams - parameters for likelihood function. In order: 
%       Sigma - variance of likelihood function)
%       PosBnds - upper and lower bounds for block position, e.g., 0 and 10 
%           for first 10 trials
%       itiBnds - upper and lower percentils for lat2fit data to be
%       included (e.g., 10 and 90 would exclude data below 10th percentile
%           and above 90th percentile)
% OUTPUTS:
%   nLL - Total negative log-likelihood of parameters
%   nll_Trial - Negative log-likelihood on each trial

% Generate model estimates
if strcmp(mdl, 'VanillaAlpha')
    lat_mdl = GenerateLatencyData_VanillaAlpha(params,...
        ratTrial, false, 'logn');
else
    error('Input valid model name')
end

% Pull likelihood parameters
Sigma = likparams(1);
PosBnds = [likparams(2) likparams(3)];

if length(Sigma) == 3
    itiBnds = [10 90];
else
    itiBnds = [likparams(4) likparams(5)];
end

% Nan out first trial of each session
trial_num = cell2mat(arrayfun(@(n) (1:n)', ratTrial.ntrials,...
    'UniformOutput', false));
lat2fit(trial_num == 1) = nan;

% Pull trials according to inclusion criteria
postvios = logical([0; ratTrial.vios(1:end-1)]);
usethese = postvios &...
    lat2fit > prctile(lat2fit, itiBnds(1)) &...
    lat2fit < prctile(lat2fit, itiBnds(2)) &...
    ~isnan(lat_mdl) &...
    ratTrial.block==1 &...
    ratTrial.BlockPosition >= PosBnds(1) &...
    ratTrial.BlockPosition <= PosBnds(2);

Lat2Fit = lat2fit(usethese); % Rat data
LatMdl = lat_mdl(usethese); % Model data

% sigma and mu parameters of a log-normal distibution to give a desired
% mean and variance
var = log(Sigma./LatMdl.^2 + 1);
mu = log(LatMdl) - (var/2);

% Evaluate nLL on each trial (up to a constant)
nll_trial = log(Lat2Fit) + 0.5*log(var) +...
    (log(Lat2Fit)-mu).^2./(2*var);

% Sum over trials
nLL = sum(nll_trial);

nLL_Trial = nan(size(lat2fit));
nLL_Trial(usethese) = nll_trial;
end