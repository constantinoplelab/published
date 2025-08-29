function [blockbins, all] = binBlocks_SS(A, block, numbins, transitionType, idx, trialType)
%select trials from different bins of a certain
%block type. cmc August 2024. Added first quartile up to first
%incongruent sss August 1, 2025

% INPUTS:
%   A = behavior data struct.
%   block = block type you want to select (1 = mixed, 2 = high, 3 = low).
%   numbins = how many evenly spaced bins to break the block into
%   transitiontype = the preceding block type, calculated by
%       diff(A.block). So -1 = high-mixed (2-1)
%   idx = logical vectro indicating which trials to include/exclude
%       (e.g., if you want to look early vs. late in training)
%   trialType = argument to specify whether to use all trials in q1 (0),
%       from the first incongruent trial (1), or up to the first incongruent
%       trial (2)

% OUTPUTS:
%   blockbins = cell(numbins,1) of trial indices in each bin.
%   all = combination of all bins, e.g., if you just want to condition on
%       previous block type

for j = 1:numbins
    blockbins{j} = [];
end

deltablk = [nan; diff(A.block)];
intothisblock = find(deltablk==transitionType); %e.g., mixed block following low block
all = [];
mintnum = 40; %if there are fewer than this number of trials in the block, don't bother trying to break it into quartiles


if transitionType==-2
    incongr = find(A.block==1 & A.reward>3); %post-low block.
elseif transitionType==-1
    incongr = find(A.block==1 & A.reward<3); %post-high block.
end

new = find(A.block~=block);
session_transition = [find(diff(A.trial_num)~=1); length(A.reward)];
    
for j = 1:length(intothisblock) %go through each block transition
    startt = intothisblock(j);
    endt1 = session_transition(min(find(session_transition>startt)));
    endt2 = new(min(find(new>startt)))-1; %just find first trial in next block

    if (endt2 - startt) < mintnum %only break into bins if there are enought trials
        continue
    end

    it = incongr(min(find(incongr >= startt)));
    if trialType==1
        startt = it+1;
    end
    
    if trialType == 2
        if it > startt
            endt2 = it - 1;
        elseif it == startt
            endt2 = nan;
        else 
            endt2 = endt1;
        end
    end
    
    if endt2>endt1 %if trials exceed the end of the session, nan it.
        endt2 = nan;
        endt1 = nan;
    end

    endt = min([endt1 endt2]);

    if ~isnan(endt1) & ~isnan(endt2) & ~isempty(startt)
        trials = startt:endt;

        %get rid of trials that are not in the block type of interest.
        [~, ia] = intersect(trials, new);
        trials(ia) = nan;
        trials(isnan(trials)) = [];

        [~, ia] = setdiff(trials, find(idx));
        if ~isempty(ia)
            trials(ia) = [];
        end

        if trialType == 2
            blockbins{1} = [blockbins{1}; (startt:endt2)']; %trials up to first incongruent
        else
             if ~isempty(trials) %& length(trials)>mintnum
                all = [all; trials'];
                %break up block into equally spaced bins
                bins = linspace(trials(1), trials(end), numbins+1);
                for kk = 2:length(bins)
                    if kk==2 %be inclusive of the first trial.
                        blockbins{kk-1} = [blockbins{kk-1}; (trials(trials>=bins(kk-1) & trials<=bins(kk)))'];
                    else
                        blockbins{kk-1} = [blockbins{kk-1}; (trials(trials>bins(kk-1) & trials<=bins(kk)))'];
                    end
                end
            end
        end
    end

end
