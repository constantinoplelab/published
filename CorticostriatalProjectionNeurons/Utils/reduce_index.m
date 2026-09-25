function [all_SU,all_S,all_index] = reduce_index(index,S, SU, goods)
% reduces the size of the SU struct to only include a subset of cells. 


all_SU = {};
all_S = {};
%all_pt = {};
%all_behEvents = {};
all_index = [];

for i = 1:length(goods)
S_goods(i) = index{goods(i),4};
end
these= unique(S_goods);

Snum = 1;
counter = 1;
for fnum = 1:length(these)
    S_temp = S{these(fnum)};
    %all_pt = [all_pt,{S_temp.pt}];
    all_S = [all_S,S_temp];
    %all_behEvents = [all_behEvents,S_temp.behEvents];


    SU_of_S = [index{find([index{:,4}]==these(fnum)),5}]';
    SU_temp = intersect(SU_of_S,goods);
    all_SU = [all_SU,SU(SU_temp)];
    
    for cluster = 1:length(SU_temp)
        all_index{counter,1} = S_temp.RatName;
        all_index{counter,2} = string(datetime(S_temp.SessionDate,'Format','yyyy-MM-dd'));
        try
            all_index{counter,3} = SU{SU_temp(cluster)}.cluster_id;
        catch
            keyboard
        end
        
        all_index{counter,4} = Snum;
        all_index{counter,5} = counter;
        counter = counter+1;
    end
    Snum = Snum+1;
end
