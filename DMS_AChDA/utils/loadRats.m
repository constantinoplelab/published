function ratlist = loadRats(datadir, sensortype)
% load list of dopamine or acetylcholine rats
% sensortype: 'da' or 'ach'

if strcmpi(sensortype, 'da')
    sensortype = 'da';
elseif strcmpi(sensortype, 'ach')
    sensortype = 'ach';
end

load(fullfile(datadir, 'data-published', 'ratlist.mat'), 'sensorList')
allrats = fields(sensorList);
ratlist = {};
for r=1:length(allrats)
    rat = allrats{r};
    if contains(sensorList.(rat), sensortype)
        ratlist = [ratlist; rat];
    end
end

end