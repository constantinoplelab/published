function processEphysData(ephysPath, codePath, savePath)
% process firing rate resoponses for 20ul trials before and after the first
% incongruent trial

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

%% load data
matfiles = dir(fullfile(ephysPath,'*.mat'));
alignto = {'COFF' 'SON' 'SOFF' 'Rew'};
ne = length(alignto);
zarg = 1;

pre_pref = cell(1, ne);
pre_nonpref = cell(1, ne);
post_pref = cell(1, ne);
post_nonpref = cell(1, ne);

info = [];

for s = 1:length(matfiles)
    load([matfiles(s).folder filesep matfiles(s).name], 'SU', 'S')
    disp(s)

    n = length(SU);

    check = find(diff(S.Block));

    if length(check) < 4 || n < 5 %require sessions to have all 3 block types and at least 5 cells
        continue
    else

        sess = string(matfiles(s).name);
        filename = repmat(sess, n, 1);

        [pre20_pref, pre20_nonpref, post20_pref, post20_nonpref, sig] = ...
            arrayfun(@(x) incongruentResponse_twentyPreVsPost(SU, S, alignto{x}, zarg), ...
            1:ne, 'UniformOutput', false);

        pre_pref = arrayfun(@(x) vertcat(pre_pref{x}, pre20_pref{x}), 1:ne, ...
            'UniformOutput', false);
        pre_nonpref = arrayfun(@(x) vertcat(pre_nonpref{x}, pre20_nonpref{x}), 1:ne, ...
            'UniformOutput', false);

        post_pref = arrayfun(@(x) vertcat(post_pref{x}, post20_pref{x}), 1:ne, ...
            'UniformOutput', false);
        post_nonpref = arrayfun(@(x) vertcat(post_nonpref{x}, post20_nonpref{x}), 1:ne, ...
            'UniformOutput', false);
        
        %keep track of which cells are significant
        id = arrayfun(@(x) SU(x).cluster_id, 1:n)';
        names = [{'id'} alignto];
        sigTable = [id cell2mat(sig)];
        i = array2table(sigTable, 'VariableNames', names);
        info = [info; addvars(i, filename, 'before', 'id')];

    end
end

% incongruentResponses.pre_pref = pre_pref;
% incongruentResponses.pre_nonpref = pre_nonpref;
% incongruentResponses.post_pref = post_pref;
% incongruentResponses.post_nonpref = post_nonpref;
% incongruentResponses.info = info;

save([savePath 'incongruentResponses.mat'], 'pre_pref', 'pre_nonpref', ...
    'post_pref', 'post_nonpref', 'info')
