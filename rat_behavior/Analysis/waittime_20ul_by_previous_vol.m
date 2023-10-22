function [wait_time_mean, pval] =...
    waittime_20ul_by_previous_vol(A, varargin)
%waittime_20ul_by_previous_vol Finds wait time on 20 uL trials in mixed 
% blocks divded into whether the previous volume is < 20 uL or > 20 uL.
%   Input: A - A struct (if you accidentally put in S struct, code will
%       convert)
%          varargin = vector of trials within the A struct that you want to
%          analyze.
%
%   Output: wait_time_mean - Average wait time for < 20 uL and > 20 uL
%       wait_time_sem - SEM for wait time for < 20 uL and > 20 uL

if isfield(A, 'pd') %just in case the input is an S struct.
    [A, ~, ~] = parse_data_from_mysql(A);
end

[~, rew] = convertreward(A.reward);
wt = A.wait_time;

if ~isempty(varargin)
    [~, ia] = setdiff(1:length(wt), varargin{1}); %set elements not in varargin to nan
    wt(ia) = nan;
end

find_20ul_catch = find(rew==20 & A.block==1 &...
    A.optout==1 & A.catch==1);
find_20ul_catch(find_20ul_catch == 1) = [];

wt = (wt - mean(wt(find_20ul_catch), 'omitnan'))./...
    std(wt(find_20ul_catch), 'omitnan');

less_than_20 = find_20ul_catch(rew(find_20ul_catch-1) < 20);
greater_than_20 = find_20ul_catch(rew(find_20ul_catch-1) > 20);

try
    pval = ranksum(wt(less_than_20), wt(greater_than_20));
catch
    pval = nan;
end

wait_time_mean = [mean(wt(less_than_20), 'omitnan'),...
    mean(wt(greater_than_20), 'omitnan')];

end

