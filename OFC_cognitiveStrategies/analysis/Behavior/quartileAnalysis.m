function [postLow, postHigh, postLow_q1, postHigh_q1] = quartileAnalysis(data)
% split mixed blocks into quartiles to look at changes in wait time over
% the course of a block

% INPUTS:
%   data = data struct containing behavior data containing wait_time,
%   reward, etc fields

data.reward = convertreward(data.reward);

%get trial indices for mixed blocks split into quartiles following low
%and high blocks.
quartiles.postLow = bin_blocks(data, 1, 4, -2, 1:length(data.wait_time), 0); %quartiles are the same for expert, inference model, and divisive normalization model
quartiles.postHigh = bin_blocks(data, 1, 4, -1, 1:length(data.wait_time), 0);

%compute wait-times by reward for the first quartile following low and high blocks
postLow_q1 = arrayfun(@(x) mean(data.wait_time(intersect(quartiles.postLow{1}, ...
    find(data.reward == x))), 'omitnan'), 1:5);
postHigh_q1 = arrayfun(@(x) mean(data.wait_time(intersect(quartiles.postHigh{1}, ...
    find(data.reward == x))), 'omitnan'), 1:5);

% compute average z-scored wait-times for each quartile in mixed blocks
% following low and high blocks
data = deltawt(data); %normalize wait times by reward to combine over rewards

postLow = arrayfun(@(x) mean(data.wait_time(intersect(quartiles.postLow{x}, ...
    find(data.optout == 1))), 'omitnan'), 1:4);
postHigh = arrayfun(@(x) mean(data.wait_time(intersect(quartiles.postHigh{x}, ...
    find(data.optout == 1))), 'omitnan'), 1:4);

end


function [A] = deltawt(A)

rew = 1:5;
for j = 1:length(rew)
    m = find(A.reward==rew(j) & A.optout==1 & A.vios==0);
    A.wait_time(m) = (A.wait_time(m)-mean(A.wait_time(m), 'omitnan'))./...
        std(A.wait_time(m), 'omitnan');
end

end