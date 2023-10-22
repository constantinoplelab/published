function [latency_mean, latency_sem, pval] =...
    latency_post20ul_by_previous_vol(A, varargin)
%latency_post20ul_by_previous_vol Finds wait time on 20 uL trials in mixed 
% blocks divded into whether the previous volume is < 20 uL or > 20 uL.
%   Input: A - A struct (if you accidentally put in S struct, code will
%       convert)
%          varargin = vector of trials within the A struct that you want to
%          analyze.
%
%   Output: latency_mean - Average wait time for < 20 uL and > 20 uL
%   cmc 11/18/21.

if isfield(A, 'pd') %just in case the input is an S struct.
    [A, ~, ~] = parse_data_from_mysql(A);
end

[~, A.reward] = convertreward(A.reward);

A.ITI(A.ITI>prctile(A.ITI, 99)) = nan;
latency = (A.ITI-mean(A.ITI, 'omitnan'))/std(A.ITI, 'omitnan'); %z-score

%set elements not in varargin to nan
if ~isempty(varargin) 
    [~, ia] = setdiff(1:length(A.ITI), varargin{1}); 
    latency(ia) = nan;
end

prevrew = [nan; A.reward(1:end-1)];
blk = A.block;

less_than_20 = prevrew<20 & blk==1;
greater_than_20 = prevrew>20 & blk==1;

latency_mean = [mean(latency(less_than_20), 'omitnan'),...
    mean(latency(greater_than_20), 'omitnan')];
latency_sem =...
    [std(latency(less_than_20), 'omitnan')./sqrt(length(less_than_20)),...
    std(latency(greater_than_20), 'omitnan')./sqrt(length(less_than_20))];

try
    pval = ranksum(latency(less_than_20), latency(greater_than_20));
catch
    pval = nan;
end


end