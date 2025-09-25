function [ms, amp] = get_avgLatency(data, type, tf)
% Get latency (ms) to first peak/trough after event onset
% 
% INPUTS
% data: 1 x 7229 event-aligned photometry data
% type (string): 'peak'/'max' or 'trough'/'min'
% tf: time window in seconds to find peak/trough. E.g., 0.5
% Default is 0.5s

T = linspace(-5, 10, 7229);

if nargin<3
    tf = 0.5;
end

[~, t0] = min(abs(T));
[~ , t1] = min(abs(T-tf));

if strcmpi(type, 'peak') || strcmpi(type, 'max')
    [amp, ms] = max(data(t0:t1)); 
else
    [amp, ms] = min(data(t0:t1));
end

% convert latency into ms
hz = 1/(T(3)-T(2));
ms = ms/hz*1000;

end