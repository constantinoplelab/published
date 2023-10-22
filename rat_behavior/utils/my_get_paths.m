function [datapath, ratList, sexList] = my_get_paths
datapath =...
    ['..' filesep 'Data' filesep 'A_Structs_Final' filesep];
a = load([datapath 'ratList.mat']);

ratList = a.ratList;
sexList = a.sexList;
end