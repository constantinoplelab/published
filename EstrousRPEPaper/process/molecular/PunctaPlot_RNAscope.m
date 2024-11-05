function TABLE = PunctaPlot_RNAscope(experiment)
%Input: experiment (eg, 'RNAscope_shEsr1validation')

%Set data path
datadir = ['Z:\histology\Carla\Confocal\' experiment ...
    '\Zstacks\compiled\Spot_data'];

%Get files to analyze
fileinfo = struct2table(dir(datadir));
fileinfo = fileinfo(3:end, :);
filenames = unique(cellfun(@(x) x, fileinfo.name, 'UniformOutput', false));
varTypes = ["string","string","string","string","double","double"];
varNames = ["RatID","Group","Region","Probe","Number_of_spots",...
    "Number_of_images"];
sz = [length(filenames) length(varNames)];
TABLE = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

for f = 1:length(filenames)
    filename = filenames{f};
    Data = readtable([datadir '\' filename]);
    underscore_idx = strfind(filename,'_');
    RatID = filename(1:underscore_idx(1)-1);    
    disp(RatID)
    Group = extractBetween(filename,'_','_');
    Region = extractBetween(filename,underscore_idx(2)+1,...
        underscore_idx(3)-1);
    Probe = extractBetween(filename,underscore_idx(end)+1, '.xlsx');
    Number_of_spots = size(Data, 1); %edited excel sheets so that all rows are spot data

    %get number of images from zstack folder
    thisfolder = ['Z:\histology\Carla\Confocal\' experiment...
        '\Zstacks\' RatID '_' Group{1} '_' Region{1}];
    thisfolder_info = struct2table(dir(thisfolder));
    uniquefilesinlist = find(contains(thisfolder_info.name, 'ch02'));
    Number_of_images = length(uniquefilesinlist);
    TABLE(f, :) = {RatID, Group, Region, Probe, Number_of_spots,...
        Number_of_images};
end

close all

end