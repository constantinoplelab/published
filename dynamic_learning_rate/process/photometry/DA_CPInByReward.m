function [CPInByRewMean, CPInByRewSEM] =...
    DA_CPInByReward(pdata, bdata)
% DA_CPInByReward - average photometry response by reward in mixed blocks
% INPUTS:
%   pdata - photometry data
%   bdata - corresponding behavioral data
% OUTPUTS:
%   CPInByRewMean - mean response
%   CPInByRewSEM - SEM response

% Preallocate data structures
[CPInByRewMean, CPInByRewSEM] = deal(nan(5, size(pdata, 2)));

% Loop over rewards
for rew = 1:5
    usethese = convertreward(bdata.Reward)==rew &...
        bdata.Block==1 &...
        bdata.PrevTrialType ~= 2;

    CPInByRewMean(rew,:) =...
        mean(pdata(usethese,:), 'omitnan');
    CPInByRewSEM(rew,:) =...
        sem(pdata(usethese,:), 'omitnan');
end
end