function [congruent_preferred, congruent_nonpreferred, ...
    incongruent_preferred, incongruent_nonpreferred, info] = ...
    incongruent_overSessions(dataPath, alignto)

%concatenate responses to incongruent responses from all cells

%INPUTS:
%   dataPath = path to .mat files containing S and SU structs for the data
%       you want to analyze
%   alignto = task events you want to align to

matfiles = dir(fullfile(dataPath,'*.mat'));

congruent_preferred = [];
congruent_nonpreferred = congruent_preferred;

incongruent_preferred = congruent_preferred;
incongruent_nonpreferred = congruent_preferred;

info = [];

for s = 1:length(matfiles)
    load([matfiles(s).folder filesep matfiles(s).name], 'SU', 'S')
    disp(s)

    n = length(SU);

    check = find(diff(S.Block) < 0);

    if length(check) < 2 || n < 5 %require sessions to have all 3 block types and at least 5 cells
        continue
    else
        sess = string(matfiles(s).name);
        filename = repmat(sess, n, 1);
        [cp, cn, ip, in, t] = ...
            incongruentRewardResponse(SU, S, alignto, 0);
        
        info = [info; addvars(t, filename, 'before', 'id')];

        congruent_preferred = vertcat(congruent_preferred, cp);
        congruent_nonpreferred = vertcat(congruent_nonpreferred, cn);
        incongruent_preferred = vertcat(incongruent_preferred, ip);
        incongruent_nonpreferred = vertcat(incongruent_nonpreferred, in);

    end
end
