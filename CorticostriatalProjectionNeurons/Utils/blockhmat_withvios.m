function [realblock, treal, tincong, vios] = blockhmat_withvios(A)
%Finds indices of block transitions for the true block transition and
%aligned to incongruent. Also outputs violation trial indices. 

tot_trials = [0; cumsum(A.ntrials)];

[~, A.reward] = convertreward(A.reward);

realblock = nan(length(A.ntrials), max(A.ntrials));
treal = realblock;
vios = realblock;
rew = realblock;
tincong = realblock;

for sess = 1:length(tot_trials)-1
    realblock(sess, 1:A.ntrials(sess)) =...
        A.block(tot_trials(sess)+1:tot_trials(sess+1));
    treal(sess,1:A.ntrials(sess)) =...
        [0; diff(A.block(tot_trials(sess)+1:tot_trials(sess+1)))];
    vios(sess,1:A.ntrials(sess)) =...
        A.vios(tot_trials(sess)+1:tot_trials(sess+1));
    rew(sess,1:A.ntrials(sess)) =...
        A.reward(tot_trials(sess)+1:tot_trials(sess+1));
end

for sess = 1:size(treal, 1)
    
    % Finds index of first mixed block trial following high block
    htom = find(treal(sess,:) == -1);
    htom_incong = nan(size(htom));
    
    for tt = 1:length(htom)
        i = find(ismember(rew(sess, htom(tt):end), [5 10]), 1)-1;
        
        if ~isempty(i)
            htom_incong(tt) = htom(tt) + i;
        else
            htom_incong(tt) = nan;
        end
    end
    htom_incong(isnan(htom_incong)) = [];
    
    % Finds index of first mixed block trial following low block
    ltom = find(treal(sess,:) == -2);
    ltom_incong = nan(size(ltom));
    
    for tt = 1:length(ltom)
        i = find(ismember(rew(sess, ltom(tt):end), [40 80]), 1)-1;
        
        if ~isempty(i)
            ltom_incong(tt) = ltom(tt) + i;
        else
            ltom_incong(tt) = nan;
        end
    end
    ltom_incong(isnan(ltom_incong)) = [];
    
    tincong(sess,1:A.ntrials(sess)) =...
        zeros(1, A.ntrials(sess));
    tincong(sess,htom_incong) = -1;
    tincong(sess,ltom_incong) = -2;
end

% Zero out block transitions before 40 trials - a bug
treal(:, 1:40) = 0;
tincong(:, 1:40) = 0;

end