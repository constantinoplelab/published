function [AUCMean, AUCSEM] =...
    DA_CPInByRPE(AUC, RPE, RPEBins, usethese, rpeBnd)
% DA_CPInByRPE - Calculates CPIn Response by RPE
% INPUTS:
%   AUC - Trial-by-trial AUC response
%   RPE - Trial-by-trial RPE
%   RPEBins - bins for RPE discretization
%   usethese - which trials to include
%   rpeBnd - upper and lower bound to exclude RPEs (symmetrical so only one
%       bound is needed)
% OUTPUTS:
%   AUCMean - Average AUC for each RPEBin
%   AUCSEM - SEM AUC for each RPEBin

% Find number of bins
nbins = length(RPEBins)-1;

% Remove RPEs above or below bound
if nargin == 4
    rpeBnd = inf;
end

RPE(abs(RPE)>=rpeBnd) = nan;

% Discritize RPEs and calculate average + sem DA response for each bin
iBin = discretize(RPE, RPEBins);
AUCMean = arrayfun(@(n)...
    mean(AUC(usethese & iBin==n), 'omitnan'), 1:nbins);
AUCSEM = arrayfun(@(n)...
    sem(AUC(usethese & iBin==n), 'omitnan'), 1:nbins);


end