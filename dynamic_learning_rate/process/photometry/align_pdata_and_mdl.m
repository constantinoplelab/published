function [RPE, AUC, Kappa] =...
    align_pdata_and_mdl(ratname, pdata, bdata,...
    BestFitEarly, BestFitLate,...
    auc1, auc2)
% align_pdata_and_mdl - aligns model fits and photometry data
% INPUTS:
%   ratname - rat ID
%   pdata - photometry data
%   bdata - corresponding behavioral data
%   BestFitEarly - Early model fit
%   BestFitLate - Late model fit
%   auc1 - time point for lower bound of AUC
%   auc2 - time point for upper bound of AUC

% Pull relevant data from bdata
ratTrial = struct('ITI', bdata.ITI, 'reward', bdata.Reward,...
    'ntrials', arrayfun(@(date) sum(bdata.UniqueDay == date),...
    unique(bdata.UniqueDay)));

% Preallocate data structures
RPE = struct;

% Get model esimated RPE for early and late
paramsEarly = BestFitEarly.(ratname).final_params(1,:);
[~, ~, RPE.Early, Kappa.Early] =...
    GenerateLatencyData_VanillaAlpha(paramsEarly, ratTrial,...
    false, 'none', nan);

paramsLate = BestFitLate.(ratname).final_params(1,:);
[~, ~, RPE.Late, Kappa.Late] =...
    GenerateLatencyData_VanillaAlpha(paramsLate, ratTrial,...
    false, 'none', nan);

% Find indices for AUC
[~, i1] = min(abs(pdata.times - auc1));
[~, i2] = min(abs(pdata.times - auc2));
dt = mean(diff(pdata.times)); % sample frequency

% Calculate AUC
AUC = sum(pdata.Data(:, i1:i2), 2)*dt;

end