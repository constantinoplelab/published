function bstruct = getBstruct(datadir, region, sensor, ratname)

basedir = fullfile(datadir, 'data-published\PhotometryData');

[region, sensor] = standardizeInputs(region, sensor);


if strcmp(sensor, 'gACh') || strcmp(sensor, 'rDA')
    fulldatadir = fullfile(basedir, strcat('GRAB_rDAgACh_',region));
else
    fulldatadir = fullfile(basedir, strcat('GRAB_', sensor, '_', region));
end

if strcmp(sensor, 'gACh') || strcmp(sensor, 'rDA')
    bfile = fullfile(fulldatadir, strcat(ratname, '_rDAgACh_bstruct.mat'));
else
    bfile = fullfile(fulldatadir, strcat(ratname, '_', sensor, '_bstruct.mat'));
end

load(bfile, 'bstruct');
disp(bfile)


end
