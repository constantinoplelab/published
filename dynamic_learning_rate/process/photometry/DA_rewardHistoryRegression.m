function betas =...
    DA_rewardHistoryRegression(pdata, bdata, nback, usethese,...
    auc1, auc2)
% DA_rewardHistoryRegression - Regresses dopamine response against previous
% rewards
% INPUTS:
%   pdata - photometry data
%   bdata - corresponding behavioral data
%   nback - number of trials back to include
%   usethese - extra constrains on which trials to include
%   auc1 - time point for lower bound of AUC
%   auc2 - time point for upper bound of AUC
% OUTPUTS:
%   betas - regression coefficients from robustfit

% Find indicies for lower and upper bounds for AUC
[~, i1] = min(abs(pdata.times - auc1));
[~, i2] = min(abs(pdata.times - auc2));
dt = mean(diff(pdata.times));

% Find each session
dates = unique(bdata.UniqueDay);

% Preallocate data strcutres
[R, AUC] = deal(cell(length(dates), 1));

% Loop over dates
for dd = 1:length(dates)
    thisday = bdata.UniqueDay == dates(dd);

    trialNos = bdata.TrialNumber(thisday);
    R{dd} = nan(max(trialNos), 1);
    R{dd}(trialNos) = bdata.Reward(thisday);

    AUC{dd} = sum(pdata.Data(:, i1:i2), 2)*dt;
end

% Log-transform rewards
R = cell2mat(R);
R = log(R);
AUC = cell2mat(AUC);

% Make design matrix
X = nan(length(R), nback+1);
for ii = 1:nback+1
    X(:,ii) = [nan(ii-1, 1); R(1:end-(ii-1))];
end
X(isnan(X(:,1)),:) = [];

% Perform regression
betas =...
    robustfit(X(usethese,:), AUC(usethese));
end