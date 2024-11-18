function [blockbins, all] = bin_blocks(A, block, numbins, transitiontype, idx, exclude_congr)
%%%select trials from different bins of a certain
%%%block type. cmc August 2024
%% INPUTS:
%% 1. A is behavior data struct.
%% 2. block is the block type you want to select (1 = mixed, 2 = high, 3 = low).
%% 3. transitiontype specifies the preceding block type, calculated by
%% diff(A.block). So -1 = high-mixed (2-1)
%% 4. idx is a logical vector indicating which trials to include/exclude
%% (e.g., if you want to look early vs. late in training)
%% 5. exclude_congr is a logical that indicates whether the mixed block
%% should be divided into terciles using all trials, or just those following
%% the first incongruent one.
%% OUTPUTS:
%% 1. blockbins is a cell(numbins,1) of trial indices in each bin.
%% 2. all combines all bins, e.g., if you just want to condition on
%% previous block type

rew = convertreward(A.reward);
for j = 1:numbins
    blockbins{j} = [];
end

deltablk = [nan; diff(A.block)];
intothisblock = find(deltablk==transitiontype); %e.g., mixed block following low block
all = [];
mintnum = 40; %if there are fewer than this number of trials in the block,
% don't bother trying to break it up into terciles.

if exclude_congr==1
    if transitiontype==-2
        incongr = find(A.block==1 & rew>3); %post-low block.
    elseif transitiontype==-1
        incongr = find(A.block==1 & rew<3); %post-high block.
    end
end

for j = 1:length(intothisblock) %go through each block transition
    startt = intothisblock(j);
    new = find(A.block~=block);
    session_transition = find(diff(A.trial_num)~=1);
    endt1 = session_transition(min(find(session_transition>startt)));
    endt2 = new(min(find(new>startt)))-1; %just find first trial in next block
    if exclude_congr==1
        startt = incongr(min(find(incongr>startt)))+1;
    end
    if endt2>endt1 %if trials exceed the end of the session, nan it.
        endt2 = nan;
        endt1 = nan;
    end

    endt = min([endt1 endt2]);
    if ~isnan(endt1) & ~isempty(startt)
        trials = startt:endt;

        %early in training there might be some trials in training stage 8
        [~, ia] = intersect(trials, find(A.trainingstage~=9));
        trials(ia) = nan;

        %get rid of trials that are not in the block type of interest.
        [~, ia] = intersect(trials, find(A.block~=block));
        trials(ia) = nan;
        trials(isnan(trials)) = [];

        [~, ia] = setdiff(trials, find(idx));
        if ~isempty(ia)
            trials(ia) = [];
        end

        if ~isempty(trials) & length(trials)>mintnum
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