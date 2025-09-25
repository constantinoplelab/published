function [pdata_avg, pdata] = get_avgFl(datadir, region, sensor, datatype)

basedir = fullfile(datadir, 'data-published\PhotometryDataAligned');
[region, sensor] = standardizeInputs(region, sensor);
fulldatadir = fullfile(basedir, strcat(sensor, '_', region));

if strcmpi(datatype, 'volume')
    namepattern = '*_avgFbyVol.mat';
elseif strcmpi(datatype, 'delay')
    namepattern = '*_avgFbyDelay.mat';
elseif strcmpi(datatype, 'side')
    namepattern = '*_avgFbySide.mat';
end

files = dir(fullfile(fulldatadir, namepattern));
files = {files.name};
fprintf('%i files found\n', length(files))

pdata = {};

for f=1:length(files)
    if ~strcmpi(datatype, 'side')
        load(fullfile(fulldatadir, files{f}), 'data');
    else
        load(fullfile(fulldatadir, files{f}), 'dataContra');
        data = dataContra;
    end

    for r=1:size(data,1)
        for c=1:size(data,2)
            pdata{r,c}(f,:) = data{r,c};
        end
    end
end 

pdata_avg = cellfun(@(x) mean(x,1,'omitnan'), pdata, UniformOutput=false);
