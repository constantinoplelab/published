function [ratBlk, mdlBlk] =...
    compare_rat_and_model(BestFit)
% compare_rat_and_model - function to compare rat trial intiation time
% and corresponding model predictions
% INPUTS:
%   BestFit - Fit to rat's data, e.g., BestFitLate.(RAT_ID)
% Outputs:
%   ratBlk - First row is average rat trial initiation time by block, 
%       second is SEM. Column order is Mixed, High, Low
%   mdlBlk - Model predicted trial intiation time. Same format as rat

% Pull unseen test data
A = BestFit.Test.ratTrial;

% Remove first trial 
trial_num = cell2mat(arrayfun(@(n) (1:n)', A.ntrials,...
    UniformOutput=false));
A.ITI(trial_num == 1) = nan;

% Only include post-violation trials (see Mah at al. 2023)
postvios = logical([0; A.vios(1:end-1)]);

% Pull 1. post-violation trials with 2. ITIs within bounds that are 3.
% Within the last 10 trials of each block
usethese = postvios &...
    A.ITI < prctile(A.ITI, 90) &...
    A.BlockPosition >= -10 &...
    A.BlockPosition <= 0;

ITIRat = A.ITI(usethese); % rat data
ITIMdl = BestFit.Test.LatMdl(usethese); % model data
Blk = A.block(usethese); % block vector

% Average trial intiation time by block
ratBlk(1,:) = arrayfun(@(b) mean(ITIRat(Blk == b), 'omitnan'), 1:3);
mdlBlk(1,:) = arrayfun(@(b) mean(ITIMdl(Blk == b), 'omitnan'), 1:3);

% SEM 
ratBlk(2,:) =...
    arrayfun(@(b)...
    std(ITIRat(Blk == b), 'omitnan')./sqrt(sum(Blk==b)), 1:3);
mdlBlk(2,:) =...
    arrayfun(@(b)...
    std(ITIMdl(Blk == b), 'omitnan')./sqrt(sum(Blk==b)), 1:3);
end