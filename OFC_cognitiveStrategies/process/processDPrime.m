function processDPrime(ephysPath, codePath, savePath)
% Compute the discriminability index for large versus small volume rewards
% in mixed blocks for each cell in the recorded area of interest

%INPUTS: 
%   ephysPath = local path to ephys folder downloaded from zenodo for the
%       recording area you want to analyze
%   codePath = local path to where code was saved
%   savePath = local path to save output

%% Set up paths
s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codePath, 'IgnoreCase', ispc);

if ~onPath % only add code dir to path if it already isn't
    addpath(genpath(codePath))
end

%% general functions
sem = @(x) std(x,'omitnan')./sqrt(size(x,1));
setDefaultFigProps
fsize = [10 10 16 10];
excludeVio = @(x, y) intersect(x, y);

alignto = {'COFF' 'SON' 'SOFF' 'Rew' 'Opt'};
ne = length(alignto);
wndw = [-0.5 4];


%%
matfiles = dir(fullfile(ephysPath, '*.mat'));

d_block = cell(1, ne);
d_vol = cell(1, ne);


for s = 1:length(matfiles)
    load([matfiles(s).folder filesep matfiles(s).name], 'SU', 'S')
    disp(s)

    n = length(SU);
    
    if n > 5
        [v, ~] = convertreward(S.RewardAmount);
        NoVio = find(S.vios == 0);
        
        % find indices of small volume (5, 10ul) and large volume (40,
        % 80ul) trials in mixed blocks)
        mix = excludeVio(find(S.Block == 1), NoVio);
        vols = arrayfun(@(x) intersect(mix, find(v == x)), ...
            1:5, 'UniformOutput', false);
    
        trials_vol{1} = [vols{1}; vols{2}]; %small rewards
        trials_vol{2} = [vols{4}; vols{5}]; %large rewards
    
        volChecks = arrayfun(@(x) length(trials_vol{x}), 1:length(trials_vol));
        
        if all(volChecks > 5)  %include sessions with at least 5 cells and at least 5 trials of each type
            [dv, ~] = arrayfun(@(x) compute_dPrime(SU, alignto{x}, wndw, ...
                    trials_vol{1}, trials_vol{2}), 1:ne, 'uniformoutput', false);
            d_vol = arrayfun(@(x) vertcat(d_vol{x}, dv{x}), 1:ne, ...
                'UniformOutput', false);
        end
    end
        
end

save([savePath 'd'], 'd_vol');
