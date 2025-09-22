function [pstruct, channel] = getPstruct(datadir, region, sensor, ratname)

basedir = fullfile(datadir, 'data-published\PhotometryData');
[region, sensor] = standardizeInputs(region, sensor);

if strcmp(sensor, 'gACh') || strcmp(sensor, 'rDA')
    fulldatadir = fullfile(basedir, strcat('GRAB_rDAgACh_',region));
else
    fulldatadir = fullfile(basedir, strcat('GRAB_', sensor, '_', region));
end

namepattern = strcat(ratname, '_', sensor, '_ch*_', region, '.mat');
pfile = dir(fullfile(fulldatadir, namepattern));
load(fullfile(fulldatadir, pfile.name), 'pstruct');
disp(fullfile(fulldatadir, pfile.name))
where = strfind(pfile.name, '_');
channel = str2double(pfile.name(where(2)+3));

end