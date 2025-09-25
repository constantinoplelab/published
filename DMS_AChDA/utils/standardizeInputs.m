function [region, sensor] = standardizeInputs(region, sensor)

if strcmpi(region, 'dms')
    region = 'DMS';
elseif strcmpi(region, 'dls')
    region = 'DLS';
elseif strcmpi(region, 'nacc')
    region = 'NAcc';
end
if strcmpi(sensor, 'da')
    sensor = 'DA';
elseif strcmpi(sensor, 'ach')
    sensor = 'ACh';
elseif strcmpi(sensor, 'gach')
    sensor = 'gACh';
elseif strcmpi(sensor, 'rda')
    sensor = 'rDA';
end

end