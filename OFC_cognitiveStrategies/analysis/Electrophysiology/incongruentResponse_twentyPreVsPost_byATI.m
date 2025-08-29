function [low, high] = incongruentResponse_twentyPreVsPost_byATI(ephysPath, tablePath, prc)
%Look at responses to 20ul trials pre and post incongruent trial split by
%ATI

%INPUTS:
%   ephysPath = path to ephys data structs
%   tablePath = path to table with blockiness metric
%   prc = percentile to look for blockiness metric defining "low
%   blockiness"

alignto = {'COFF' 'SON' 'SOFF' 'Rew'};
nAlign = length(alignto);
zarg = 1;

matfiles = dir(fullfile(ephysPath,'*.mat'));

T = readtable(tablePath);

ATI = T.ATI;
group = string(T.Group);
ATI = ATI(strcmp(group, 'Naive')); %only use naive rat sessions

lowATI = find(ATI < prctile(ATI, prc));
highATI = find(ATI > prctile(ATI, 100-prc));

sessions = T.FileName(strcmp(group, 'Naive'));

prePref_low = cell(1, nAlign);
preNonpref_low = prePref_low;
postPref_low = prePref_low;
postNonpref_low = prePref_low;

prePref_high = cell(1, nAlign);
preNonpref_high = prePref_high;
postPref_high = prePref_high;
postNonpref_high = prePref_high;

info_low = [];
for ii = 1:length(lowATI)
    
    isSess = find(arrayfun(@(x) contains(matfiles(x).name, ...
        sessions(lowATI(ii))), 1:length(matfiles)));

    load([matfiles(isSess).folder filesep matfiles(isSess).name], 'SU', 'S')
    n = length(SU);
    sess = string(matfiles(isSess).name);
    filename = repmat(sess, n, 1);
    
    [pre20_pref, pre20_nonpref, post20_pref, post20_nonpref, sig] = ...
        arrayfun(@(x) incongruentResponse_twentyPreVsPost(SU, S, alignto{x}, ...
        zarg), 1:nAlign, 'uniformoutput', false);

    prePref_low = arrayfun(@(x) vertcat(prePref_low{x}, pre20_pref{x}), 1:nAlign, ...
        'UniformOutput', false);
    preNonpref_low = arrayfun(@(x) vertcat(preNonpref_low{x}, pre20_nonpref{x}), 1:nAlign, ...
        'UniformOutput', false);

    postPref_low = arrayfun(@(x) vertcat(postPref_low{x}, post20_pref{x}), 1:nAlign, ...
        'UniformOutput', false);
    postNonpref_low = arrayfun(@(x) vertcat(postNonpref_low{x}, post20_nonpref{x}), 1:nAlign, ...
        'UniformOutput', false);

    %keep track of which cells are significant
    id = arrayfun(@(x) SU(x).cluster_id, 1:n)';
    names = [{'id'} alignto];
    sigTable = [id cell2mat(sig)];
    i = array2table(sigTable, 'VariableNames', names);
    info_low = [info_low; addvars(i, filename, 'before', 'id')];

    
end

info_high = [];
for jj = 1:length(highATI)
    
    isSess = find(arrayfun(@(x) contains(matfiles(x).name, ...
        sessions(highATI(jj))), 1:length(matfiles)));

    load([matfiles(isSess).folder filesep matfiles(isSess).name], 'SU', 'S')
    n = length(SU);
    sess = string(matfiles(isSess).name);
    filename = repmat(sess, n, 1);
    
    [pre20_pref, pre20_nonpref, post20_pref, post20_nonpref, sig] = ...
        arrayfun(@(x) incongruentResponse_twentyPreVsPost(SU, S, alignto{x}, ...
        zarg), 1:nAlign, 'uniformoutput', false);

    prePref_high = arrayfun(@(x) vertcat(prePref_high{x}, pre20_pref{x}), 1:nAlign, ...
        'UniformOutput', false);
    preNonpref_high = arrayfun(@(x) vertcat(preNonpref_high{x}, pre20_nonpref{x}), 1:nAlign, ...
        'UniformOutput', false);

    postPref_high = arrayfun(@(x) vertcat(postPref_high{x}, post20_pref{x}), 1:nAlign, ...
        'UniformOutput', false);
    postNonpref_high = arrayfun(@(x) vertcat(postNonpref_high{x}, post20_nonpref{x}), 1:nAlign, ...
        'UniformOutput', false);
    
    %keep track of which cells are significant
    id = arrayfun(@(x) SU(x).cluster_id, 1:n)';
    names = [{'id'} alignto];
    sigTable = [id cell2mat(sig)];
    i = array2table(sigTable, 'VariableNames', names);
    info_high = [info_high; addvars(i, filename, 'before', 'id')];

end

low.pre_pref = prePref_low;
low.pre_nonpref = preNonpref_low;
low.post_pref = postPref_low;
low.post_nonpref = postNonpref_low;
low.info = info_low;

high.pre_pref = prePref_high;
high.pre_nonpref = preNonpref_high;
high.post_pref = postPref_high;
high.post_nonpref = postNonpref_high;
high.info = info_high;
