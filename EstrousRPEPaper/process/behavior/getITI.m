function [iti,iti_npv,iti_pv,startTime] = getITI(S)
% Calculates inter-trial interval and start time for each trial in a
% session. SSS 05/2020
%INPUTS: 
%   S = screened sessions from get_data (using poolRats to screen for 
%       current criteria 05/2020)
%OUTPUTS:
%   iti = cell for each session containing the ITI for each trial. ITI is
%       defined as side out to center in (if the previous trial was
%       rewarded and licking bleeds into the current trial, that time is
%       subtracted from CenterIn in S.peh)
%   iti_npv = cell for each session. iti with post violation trials nan'd
%   startTime = cell for each session. Time of trial start relative to
%       session start


% itiBinned = {};
iti = cell(1,length(S.pd));
iti_npv = iti; %exclude post violation trials
iti_pv = iti; %include only post violation trials
startTime = iti;
    
for s = 1:length(S.pd)

    numTrials = length(S.peh{s});
    I = NaN(numTrials+1,1);
    addTime = zeros(numTrials+1,1);

    exclude = find(S.pd{s}.vios) + 1; %exclude post violation ITIs due to vio TO

    R = intersect(find(strcmp(S.pd{s}.RewardedSide(1:numTrials),'R')), ...
        find(S.pd{s}.hits));
    L = intersect(find(strcmp(S.pd{s}.RewardedSide(1:numTrials),'L')), ...
        find(S.pd{s}.hits));
    
    if isfield(S.peh{s}(1).Events,'CenterIn')
        I(1) = S.peh{s}(1).Events.CenterIn(1);
    end
    
    %check for pokes on previous reward port before center in
    %add licking time to trial start time for current trial

    for t = 2:numTrials
        if isfield(S.peh{s}(t).Events, 'Tup')
            if ~isfield(S.peh{s}(t).Events, 'CenterIn')
                I(t) = NaN;
            else
                flag = 0;

                cIn = S.peh{s}(t).Events.CenterIn(1);
                Rlicks = [];
                Llicks = [];

                if isfield(S.peh{s}(t).Events,'RightOut')
                    rOut = S.peh{s}(t).Events.RightOut;
                    Rlicks = rOut(rOut < cIn);
                end
                if isfield(S.peh{s}(t).Events,'LeftOut')
                    lOut = S.peh{s}(t).Events.LeftOut;
                    Llicks = lOut(lOut < cIn);
                end

                if ismember(t-1,R) && ~isempty(Rlicks)
                    Rlicks = Rlicks(end); %rat finishes collecting reward from prev trial
                    I(t) = cIn - Rlicks;
                    addTime(t) = Rlicks; %actual start of trial
                elseif ismember(t-1,L) && ~isempty(Llicks)
                    Llicks = Llicks(end); 
                    I(t) = cIn - Llicks;  
                    addTime(t) = Llicks;
                else
                    I(t) = cIn;
                end
            end
        else 
            flag = 1;
            break 
        end
    end

    npv = I; npv(exclude) = nan; pv = I;  pv(~exclude) = nan;
    npv(end) = [];
    I(end) = [];
    pv(end) = [];

    iti{s} = I;
    iti_npv{s} = npv;
    iti_pv{s} = pv;


end
