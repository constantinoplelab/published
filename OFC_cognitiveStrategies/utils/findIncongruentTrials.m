function [t_incongruent, t_congruent] = findIncongruentTrials(S)
%Finds the index for the first incongruent and congruent rewarded trials

% find incongruent trials
bchange = diff(S.Block);
allChange = find(bchange);
ltom = find(bchange == -2)+1; %first mixed trial post low
ltom_incong_rew = nan(size(ltom));
ltom_cong_rew = nan(size(ltom));

htom = find(bchange == -1)+1; %first mixed trial post high
htom_incong_rew = nan(size(htom));
htom_cong_rew = nan(size(htom));

v = convertreward(S.RewardAmount);
lowRew = find(v==1 | v==2);
highRew = find(v==4 | v==5);
twenty = find(v==3);

for tt = 1:length(htom)
    next = [allChange(allChange > htom(tt)); length(v)];
    
    i = intersect(lowRew, htom(tt):next(1));
    rewarded = intersect(i, find(S.hits));
    
    try
        c = intersect([highRew; twenty], htom(tt):i(1));
    catch
        c = htom(tt):next(1); %no incongruent trials
    end
    c_rew = intersect(c, find(S.hits));
        
    if ~isempty(rewarded)
        htom_incong_rew(tt) = rewarded(1);
    else
        htom_incong_rew(tt) = nan;
    end

    if ~isempty(c_rew)
        htom_cong_rew(tt) = c_rew(1);
    else
        htom_cong_rew(tt) = nan;
    end
end

for tt = 1:length(ltom)
    next = [allChange(allChange > ltom(tt)); length(v)];

    i = intersect(highRew, ltom(tt):next(1));
    rewarded = intersect(i, find(S.hits));
    opt = intersect(i, find(S.optout));
    
    try
        c = intersect([lowRew; twenty], ltom(tt):i(1));
    catch
        c = ltom(tt):next(1); %if no incongruent trials
    end
    c_rew = intersect(c, find(S.hits));

    if ~isempty(rewarded)
        ltom_incong_rew(tt) = rewarded(1);
    else
        ltom_incong_rew(tt) = nan;
    end

    if ~isempty(c_rew)
        ltom_cong_rew(tt) = c_rew(1);
    else
        ltom_cong_rew(tt) = nan;
    end

end


ltom_incong_rew(isnan(ltom_incong_rew)) = [];
htom_incong_rew(isnan(htom_incong_rew)) = [];

t_incongruent.ltom = [ltom_incong_rew v(ltom_incong_rew)];
t_incongruent.htom = [htom_incong_rew v(htom_incong_rew)];

ltom_cong_rew(isnan(ltom_cong_rew)) = [];
htom_cong_rew(isnan(htom_cong_rew)) = [];

t_congruent.ltom = [ltom_cong_rew v(ltom_cong_rew)];
t_congruent.htom = [htom_cong_rew v(htom_cong_rew)];
