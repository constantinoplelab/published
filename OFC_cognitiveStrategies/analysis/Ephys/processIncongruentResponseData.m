function processIncongruentResponseData(dataPath, savePath)
% Pull out and save responses to congruent and incongruent trials for all
% sessions. Naming conventions used in the function refer to the block
% preference, not the transition type as used in the paper (e.g.
% congruent_preferred refer to responses to congruent trials after the
% preferred block, a non-preferred transition in the paper)

%INPUTS: 
%   dataPath = path to the ephys data you want to process (e.g. /Ephys/Expert) 
%   savePath = where you want to save the data. please save each data
%       set in its own folder (e.g. create separate folders to save expert
%       and naive data in)

%OUTPUTS
%   saves a .mat file containing the following variables
%       incongPref = responses to incongruent trials following a transition
%           from the cell's preferred block into a mixed block. Cells that 
%           do not pass block significance criterion are nan'd
%       incongNonpref = same as incongPref but for responses to transitions
%           from the cell's non-preferred block into a mixed block
%       congPref = same as incongPref but for responses to congruent trials
%       congNon = same as incongNonpref bur for responses to congruent
%           trials
     
alignto = {'COFF' 'SON' 'SOFF' 'Rew'};
nAlign = length(alignto);

[congruent_preferred, congruent_nonpreferred, ...
    incongruent_preferred, incongruent_nonpreferred, info] = ...
    incongruent_overSessions(dataPath, alignto);

incongPref = arrayfun(@(x) cell2mat(incongruent_preferred(:,x)), ...
    1:nAlign, 'UniformOutput', false);

incongNonpref = arrayfun(@(x) cell2mat(incongruent_nonpreferred(:,x)), ...
    1:nAlign, 'uniformoutput', false);

congPref = arrayfun(@(x) cell2mat(congruent_preferred(:,x)), ...
    1:nAlign, 'uniformoutput', false);

congNonpref = arrayfun(@(x) cell2mat(congruent_nonpreferred(:,x)), ...
    1:nAlign, 'uniformoutput', false);

save([savePath 'incongruentResponses.mat'], 'incongPref', 'incongNonpref', ...
    'congPref', 'congNonpref', 'info')
