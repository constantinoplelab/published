function [CPIn20ByBlkMean, CPIn20ByBlkSEM] = DA_CPInByBlock(pdata, bdata)
% DA_CPInByBlock - average photometry response by 20 ul across blocks
% INPUTS:
%   pdata - photometry data
%   bdata - corresponding behavioral data
% OUTPUTS:
%   CPIn20ByBlkMean - mean response
%   CPIn20ByBlkSEM - SEM response

% Preallocate data strucutre
[CPIn20ByBlkMean, CPIn20ByBlkSEM] = deal(nan(3, size(pdata, 2)));

% Loop over blocks
for blk = 1:3
    usethese = convertreward(bdata.Reward)==3 &...
        bdata.Block==blk &...
        bdata.PrevTrialType ~= 2;

    CPIn20ByBlkMean(blk,:) =...
        mean(pdata(usethese,:), 'omitnan');
    CPIn20ByBlkSEM(blk,:) =...
        sem(pdata(usethese,:), 'omitnan');
end

end
