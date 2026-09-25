function clusternames = getclusternames(SU,varargin)

if isempty(varargin)
for i = 1:length(SU)
    clusternames(i,1) = SU{1,i}.cluster_id;
end
else
    goods = varargin{1};
    for i = 1:length(varargin{1})
         clusternames(i) = SU{goods(i)}.cluster_id;
    end
end

end
