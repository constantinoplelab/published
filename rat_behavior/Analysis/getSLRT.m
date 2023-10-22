function [slrt] = getSLRT(S)
% Calculates side led reaction time (side led on ~ side poke) 
% Hs 2022/11
%
%INPUTS: 
%   S = screened sessions from get_data saved in
%   server/ProcessedRatData/S_Structs folder. 
%OUTPUTS:
%   slrt = cell for each session containing the side led reaction time for each trial.

slrt = cell(1,length(S.pd));

for s = 1:length(S.pd)
    numTrials = length(S.peh{s});
    I = NaN(numTrials,1);

    for t = 1:numTrials
        if isfield(S.peh{s}(t).Events, 'CenterIn') && isfield(S.peh{s}(t).Events, 'CenterOut')
            if ~isnan(S.peh{s}(t).States.WaitForSidePoke(1)) % hit or opt-out trials: get the first side poke from S.peh{s}(t).Events
                if strcmp(S.pd{s}.RewardedSide{t}, 'R')
                    try
                        xx = find(S.peh{s}(t).Events.RightIn > S.peh{s}(t).States.WaitForSidePoke(1),1); % get the first rpoke after ROn
                        I(t) = S.peh{s}(t).Events.RightIn(xx) - S.peh{s}(t).States.WaitForSidePoke(1);
                    catch
                        I(t) = NaN; % e.g., if rat opts out directly without checking out reward port
                    end
                elseif strcmp(S.pd{s}.RewardedSide{t}, 'L')
                    try
                        xx = find(S.peh{s}(t).Events.LeftIn > S.peh{s}(t).States.WaitForSidePoke(1),1); % get the first lpoke after LOn
                        I(t) = S.peh{s}(t).Events.LeftIn(xx) - S.peh{s}(t).States.WaitForSidePoke(1);
                    catch
                        I(t) = NaN; % e.g., if rat opts out directly without checking out reward port
                    end
                end
            else % violation trials
                I(t) = NaN;
            end
        else
            I(t) = NaN;
        end
    end

    slrt{s} = I;

end


end