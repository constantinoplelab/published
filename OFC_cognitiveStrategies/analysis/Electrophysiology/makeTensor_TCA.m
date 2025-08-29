function makeTensor_TCA(dataPath, savePath)
%INPUTS:
%   dataPath = path to matfiles containing S and SU structs -- these should
%       be files with only cells from the specific region of interest
%   savePath = path to folder for saving 3d matrix of responses

%OUTPUTS:
%   data = 3d matrix that will be used to create tensor (neurons X time X block)
%   tensorInfo = table with information about each cell (filename, cell Id)

matfiles = dir(fullfile(dataPath,'*.mat'));
events = {'COFF' 'SON' 'SOFF' 'Rew' 'Opt'};

nsess = length(matfiles);

data = [];
tensorInfo = [];

excludeVio = @(x, y) intersect(x, y);

for s = 1:nsess
    disp(s)
    
    name = matfiles(s).name;
    load([matfiles(s).folder filesep name], 'SU', 'S')

    n = length(SU);
    
    check = find(diff(S.Block)); %all 3 full blocks

    if n < 5 || length(check) < 4
        continue
    end
    
    NoVio = find(S.vios == 0);

    blockVec = [3, 1, 2]; %put in low, mix, high order for tensor construction
    all = arrayfun(@(x) excludeVio(find(S.Block == x), NoVio), ...
        blockVec, 'uniformoutput', false);

    data = cat(1, data, sessTensor(all, SU, events));

    sess = string(name);
    filename = repmat(sess, n, 1);
    id = arrayfun(@(x) SU(x).cluster_id, 1:n)';
    tensorInfo = [tensorInfo; table(filename, id)];

end

tensorInfo.Properties.VariableNames = ["filename", "id"];
save(strcat(savePath,'time', 'block', 'Tensor'),'data')
save(strcat(savePath, 'block', '_tensorInfo'),'tensorInfo')

end

function out = sessTensor(trials, SU, events)
%average over trial types to create session tensor
%INPUTS:
%   trials = cell array of trial indices for the different conditions to
%       average over in the z-dimension of the tensor (e.g. blocks)

zscore = @(x) (x - mean(x,'all','omitnan')) / std(x,0,'all','omitnan');

timeBin = [-1 1];
w = arrayfun(@(x) find(round(SU(1).xvec.COFF, 3) == x), timeBin);



out = [];
for e = 1:length(events)
    if strcmpi(events{e}, 'Rew') || strcmpi(events{e}, 'Opt')
        timeBin = [-1 0.8];
        w = arrayfun(@(x) find(round(SU(1).xvec.COFF, 3) == x), timeBin);
    end
    zSU = arrayfun(@(x) zscore(SU(x).hmat.(events{e})), 1:length(SU),...
        'uniformoutput', false);
    initial = zeros(length(zSU), length(w(1):w(2)), 3);
    
    for ii = 1:length(trials)
       avg_resp = arrayfun(@(x) mean(zSU{x}(trials{ii},w(1):w(2)), 'omitnan'),...
           1:length(zSU), 'uniformoutput',false); 
       avg_resp = cell2mat(arrayfun(@(x) avg_resp{x}(:), 1:length(avg_resp), ...
           'uniformoutput', false))'; 
       if ~isnan(avg_resp)
           initial(:,:,ii) = avg_resp;
       end
    end
    
    out = [out initial];
end

end