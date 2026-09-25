function [postLow, postHigh, postLow_q1, postHigh_q1, postLow_q1Con, ...
    postHigh_q1Con] = quartileAnalysis_SS(data)
% split mixed blocks into quartiles to look at changes in wait time over
% the course of a block

% INPUTS:
%   data = data struct containing behavior data containing wait_time,
%   reward, etc fields

check = find(data.reward == 1);

if isempty(check)
    data.reward = convertreward(data.reward); %convert to 1:5 if in standard format
end

%get trial indices for mixed blocks split into quartiles following low
%and high blocks.
quartiles.postLow = binBlocks_SS(data, 1, 4, -2, 1:length(data.wait_time), 1); %from first incongruent trial
quartiles.postHigh = binBlocks_SS(data, 1, 4, -1, 1:length(data.wait_time), 1);

%compute wait-times by reward for the first quartile following low and high blocks
use = arrayfun(@(x) find(data.reward == x & data.catch & data.optout ...
    & data.block == 1), ...
    1:5, 'UniformOutput', false);
postLow_q1 = arrayfun(@(x) mean(data.wait_time(intersect(quartiles.postLow{1}, ...
    use{x})), 'omitnan'), 1:5);
postHigh_q1 = arrayfun(@(x) mean(data.wait_time(intersect(quartiles.postHigh{1}, ...
    use{x})), 'omitnan'), 1:5);

%compute average wait time for q1 before the first incongruent trial
quartiles.postLow = binBlocks_SS(data, 1, 4, -2, 1:length(data.wait_time), 2); %up to first incongruent
quartiles.postHigh = binBlocks_SS(data, 1, 4, -1, 1:length(data.wait_time), 2);

postLow_q1Con = arrayfun(@(x) mean(data.wait_time(intersect(quartiles.postLow{1}, ...
    use{x})), 'omitnan'), 1:5);
postHigh_q1Con = arrayfun(@(x) mean(data.wait_time(intersect(quartiles.postHigh{1}, ...
    use{x})), 'omitnan'), 1:5);

% compute average z-scored wait-times for each quartile in mixed blocks
% following low and high blocks
quartiles.postLow = binBlocks_SS(data, 1, 4, -2, 1:length(data.wait_time), 0); %include congruent trials
quartiles.postHigh = binBlocks_SS(data, 1, 4, -1, 1:length(data.wait_time), 0);

% data = detrendwt_SS(data);
data = deltawt(data); %normalize wait times by reward to combine over rewards

postLow = arrayfun(@(x) mean(data.wait_time(intersect(quartiles.postLow{x}, ...
    find(data.optout & data.catch))), 'omitnan'), 1:4);
postHigh = arrayfun(@(x) mean(data.wait_time(intersect(quartiles.postHigh{x}, ...
    find(data.optout & data.catch))), 'omitnan'), 1:4);

end


function [A] = deltawt(A)

rew = 1:5;
for j = 1:length(rew)
    m = find(A.reward==rew(j) & A.catch & A.vios==0);
    A.wait_time(m) = (A.wait_time(m)-mean(A.wait_time(m), 'omitnan'))./...
        std(A.wait_time(m), 'omitnan');
end

end