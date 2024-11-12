function avg_resp = averageResponse(dataPath, tensorInfo, timeBin, zarg)

%INPUTS:
%   dataPath: path to ephys SU structs
%   tensor info: table containing the filename, cellID, and TCA cluster for
%       each cell
%   timeBin: bin of time over which to compute average responses
%   zarg: option to z-score or use raw firing rates

%OUTPUTS:
%   avg_resp: average response for all cells 

matfiles = dir(fullfile(dataPath,'*.mat'));
sess = string(cell2mat({matfiles.name}'));

use = ismember(sess, tensorInfo.filename); %only loop through sessions that have passed criteria for inclusion in the tensor (all 3 block types, at least 5 cells)
matfiles = matfiles(use, :);

w = arrayfun(@(x) find(-4:0.05:8 == x), timeBin); %saved hmats are [-4 8s] window around each task event

zscore = @(x) (x - mean(x(:),'omitnan')) / std(x,0,'all','omitnan');

alignto = {'COFF' 'SON' 'SOFF' 'Rew' 'Opt'};
ne = length(alignto);

nc = max(tensorInfo.TCAclust);

avg_resp = cell(1, ne);

for s = 1:length(matfiles)

    name = matfiles(s).name;
    load([matfiles(s).folder filesep name], 'SU')

    good = tensorInfo.id(ismember(tensorInfo.filename, name));

    ids = extractfield(SU, 'cluster_id');
    use = ismember(ids', good);

    SU = SU(use);
    n = length(SU);

    for ii = 1:n

        if zarg
            zSU = arrayfun(@(x) zscore(SU(ii).hmat.(alignto{x})), 1:ne,...
                'uniformoutput', false);
        else
            zSU = arrayfun(@(x) SU(ii).hmat.(alignto{x}), 1:ne,...
                'uniformoutput', false); %try with raw fr
        end

        avg = arrayfun(@(x) mean(zSU{x}(:,w(1):w(2)), 'omitnan'),...
            1:ne, 'uniformoutput',false);
    
        avg_resp = arrayfun(@(x) [avg_resp{x}; avg{x}], 1:ne, ...
            'uniformoutput', false);
    end
end