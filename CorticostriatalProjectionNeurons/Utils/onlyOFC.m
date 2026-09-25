function [OFC,notOFC] = onlyOFC(all_S,cells,all_SU,all_index)

for i = 1:length(cells)
    S = all_S{all_index{cells(i),4}};
    SU = all_SU{cells(i)};
    try
    depth(i,:) = [SU.AP,SU.ML,SU.DV];
    location{i} = SU.location;
    catch
    SU = channelLocation_MD(S,{SU});
    depth(i,:) = [SU{1}.AP,SU{1}.ML,SU{1}.DV];
    location{i} = SU{1}.location;
    end
end

location = location';


OFC = cells(strcmp(location,'OFC'));
notOFC = cells(~strcmp(location,'OFC'));