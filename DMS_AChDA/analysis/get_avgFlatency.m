function [ms, amp, ratlist] = ...
    get_avgFlatency(datadir, region, sensor, datatype, minmaxtype, tf)
% Get latency (ms) to first peak/trough after event onset

if nargin<5
    tf = 0.5;
end
basedir = fullfile(datadir, 'data-published\PhotometryDataAligned');
[region, sensor] = standardizeInputs(region, sensor);
datadir = fullfile(basedir, strcat(sensor, '_', region));

if strcmpi(datatype, 'volume')
    namepattern = '*_avgFbyVol.mat';
elseif strcmpi(datatype, 'delay')
    namepattern = '*_avgFbyDelay.mat';
elseif strcmpi(datatype, 'side')
    namepattern = '*_avgFbySide.mat';
end

files = dir(fullfile(datadir, namepattern));
files = {files.name};
fprintf('%i files found\n', length(files))
ratlist = cellfun(@(x) extractBefore(x, '_'), ...
    files, 'UniformOutput', false);
ratlist = ratlist';
ms = {};
amp = {};
for f=1:length(files)
    if ~strcmpi(datatype, 'side')
        load(fullfile(datadir, files{f}), 'data');
    else
        load(fullfile(datadir, files{f}), 'dataContra');
        data = dataContra;
    end
    
    for r=1:size(data,1)
        for c=1:size(data,2)
            [ms{r,c}(f,1), amp{r,c}(f,1)] = ...
                getFlatency(data{r,c}, minmaxtype, tf);
        end
    end

end

end

function [ms, amp] = getFlatency(data, type, tf)
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

