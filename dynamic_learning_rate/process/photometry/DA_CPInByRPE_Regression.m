function [betasPos, betasNeg, rpe_by_sess, auc_by_sess,...
statsPos, statsNeg] =...
DA_CPInByRPE_Regression(bdataMdl, AUC, RPE, usethese, rpe_bnd)
% DA_CPInByRPE_Regression - Regress DA AUC against model RPE. Regresses
% separately for positive and negative RPEs in each session
% INPUTS:
%   bdataMdl - Behavorial data (same as bdata)
%   AUC - Trial-by-trial AUC response
%   RPE - Trial-by-trial RPE
%   usethese - which trials to include
%   rpeBnd - upper and lower bound to exclude RPEs (symmetrical so only one
%       bound is needed)
% OUTPUTS:
%   betasPos - Regression coefficents for positive RPE regression across
%       sessions
%   betasNeg - Regression coefficents for negative RPE regression across
%       sessions
%   rpe_by_sess - Each session's RPE
%   auc_by_sess - Each session's AUC
%   statsPos - statistics for positive RPE regression (from robustfit)
%   statsNeg - statistics for negative RPE regression (from robustfit)

% Pull each date
dates = unique(bdataMdl.UniqueDay);
[betasPos, betasNeg] =...
    deal(nan(length(dates), 2));

% Preallocate data structures
[rpe_by_sess, auc_by_sess,...
    statsPos, statsNeg] = deal(cell(length(dates), 1));

% Loop over dates
for dd = 1:length(dates)

    % Pull data for each date
    thisday = bdataMdl.UniqueDay == dates(dd);

    rpe = RPE(usethese & thisday);
    auc = AUC(usethese & thisday);

    in_rpe_bnd = abs(rpe) < rpe_bnd;

    % Positive regression - if not enough data for robustfit, skip date
    try
        [betasPos(dd,:), statsPos{dd}] =...
            robustfit(rpe(rpe > 0 & in_rpe_bnd),...
            auc(rpe > 0 & in_rpe_bnd));
    catch
    end

    % Negative regression - if not enough data for robustfit, skip date
    try
        [betasNeg(dd,:), statsNeg{dd}] =...
            robustfit(rpe(rpe < 0 & in_rpe_bnd),...
            auc(rpe < 0 & in_rpe_bnd));
    catch
    end
    
    rpe_by_sess{dd} = rpe(in_rpe_bnd);
    auc_by_sess{dd} = auc(in_rpe_bnd);
end

end