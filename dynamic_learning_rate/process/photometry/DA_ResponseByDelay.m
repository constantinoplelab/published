function [sideOnByDelay, AUC, baseline] =...
    DA_ResponseByDelay(pdata, bdata, bins,...
    auc1, auc2, bl1, bl2)
% DA_ResponseByDelay - calculates average DA response by reward delay
% INPUTS:
%   pdata - photometry data
%   bdata - corresponding behavioral data
%   bins - bins for delays
%   auc1 - time point for lower bound of AUC
%   auc2 - time point for upper bound of AUC
%   bl1 - time point for lower bound of baseline
%   bl2 - time point for upper bound of baseline
% OUTPUTS:
%   sideOnByDelay - Average DA response for each delay bin
%   AUC - Baseline corrected AUC for each delay bin
%   baseline - average baseline


    dt = mean(diff(pdata.times)); % Sample rate

    % Find indicies for lower and upper bounds for AUC
    [~, i1] = min(abs(pdata.times - auc1));
    [~, i2] = min(abs(pdata.times - auc2));

    % Find indicies for lower and upper bounds for baseline
    [~, i3] = min(abs(pdata.times - bl1));
    [~, i4] = min(abs(pdata.times - bl2));

    % mixed block rewarded trials
    these = bdata.Block == 1 & ~bdata.Catch & ~bdata.OptOut;

    % define bins based on delay distribution
    delays = bdata.RewardDelay;
    delays(delays >= 100) = nan;

    % Discretize delays
    binI = discretize(delays, bins);

    % preallocate data structures
    sideOnByDelay = nan(length(bins)-1, size(pdata.Data, 2));
    [AUC, baseline] = deal(nan(length(bins)-1, 1));

    % Loop over bins
    for i = 1:length(bins)-1
        usethese = these & binI == i;
        sideOnByDelay(i,:) = mean(pdata.Data(usethese,:), 'omitnan');

        % Calculate baseline
        baseline(i) = mean(sideOnByDelay(i,i3:i4));

        % Calculate baseline corrected AUC
        AUC(i) = sum(sideOnByDelay(i,i1:i2))*dt - baseline(i);
    end
end