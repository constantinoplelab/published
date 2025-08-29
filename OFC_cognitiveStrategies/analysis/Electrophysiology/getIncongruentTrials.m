function [incong_offer, incong_rew, incong_opt, cong_offer, ...
    cong_rew, cong_opt] = getIncongruentTrials(S)

% find incongruent trials
bchange = diff(S.Block);
allChange = find(bchange);
ltom = find(bchange == -2)+1; %first mixed trial post low
ltom_incong_offer = nan(size(ltom));
ltom_incong_rew = nan(size(ltom));
ltom_incong_opt = nan(size(ltom));
ltom_cong_offer = nan(size(ltom));
ltom_cong_rew = nan(size(ltom));
ltom_cong_opt = nan(size(ltom));

htom = find(bchange == -1)+1; %first mixed trial post high
htom_incong_offer = nan(size(htom));
htom_incong_rew = nan(size(htom));
htom_incong_opt = nan(size(htom));
htom_cong_offer = nan(size(htom));
htom_cong_rew = nan(size(htom));
htom_cong_opt = nan(size(htom));


v = convertreward(S.RewardAmount);
lowRew = find(v==1 | v==2);
highRew = find(v==4 | v==5);
twenty = find(v==3);

for tt = 1:length(htom)
    next = [allChange(allChange > htom(tt)); length(v)];
    
    i = intersect(lowRew, htom(tt):next(1));
    rewarded = intersect(i, find(S.hits));
    opt = intersect(i, find(S.optout));
    
    try
        c = intersect([highRew; twenty], htom(tt):i(1));
    catch
        c = htom(tt):next(1); %no incongruent trials
    end
    c_rew = intersect(c, find(S.hits));
    c_opt = intersect(c, find(S.optout));
        
    if ~isempty(i)
        htom_incong_offer(tt) = i(1);
    else
        htom_incong_offer(tt) = nan;
    end

    if ~isempty(rewarded)
        htom_incong_rew(tt) = rewarded(1);
    else
        htom_incong_rew(tt) = nan;
    end

    if ~isempty(opt)
        htom_incong_opt(tt) = opt(1);
    else
        htom_incong_opt(tt) = nan;
    end

    if ~isempty(c)
        htom_cong_offer(tt) = c(1);
    else
        htom_cong_offer(tt) = nan;
    end

    if ~isempty(c_rew)
        htom_cong_rew(tt) = c_rew(1);
    else
        htom_cong_rew(tt) = nan;
    end

    if ~isempty(c_opt)
        htom_cong_opt(tt) = c_opt(1);
    else
        htom_cong_opt(tt) = nan;
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
    c_opt = intersect(c, find(S.optout));

    if ~isempty(i)
        ltom_incong_offer(tt) = i(1);
    else
        ltom_incong_offer(tt) = nan;
    end

    if ~isempty(rewarded)
        ltom_incong_rew(tt) = rewarded(1);
    else
        ltom_incong_rew(tt) = nan;
    end

    if ~isempty(opt)
        ltom_incong_opt(tt) = opt(1);
    else
        ltom_incong_opt(tt) = nan;
    end

    if ~isempty(c)
        ltom_cong_offer(tt) = c(1);
    else
        ltom_cong_offer(tt) = nan;
    end

    if ~isempty(c_rew)
        ltom_cong_rew(tt) = c_rew(1);
    else
        ltom_cong_rew(tt) = nan;
    end

    if ~isempty(c_opt)
        ltom_cong_opt(tt) = c_opt(1);
    else
        ltom_cong_opt(tt) = nan;
    end
end


ltom_incong_offer(isnan(ltom_incong_offer)) = [];
htom_incong_offer(isnan(htom_incong_offer)) = [];

incong_offer.ltom = [ltom_incong_offer v(ltom_incong_offer)];
incong_offer.htom = [htom_incong_offer v(htom_incong_offer)];

ltom_incong_rew(isnan(ltom_incong_rew)) = [];
htom_incong_rew(isnan(htom_incong_rew)) = [];

incong_rew.ltom = [ltom_incong_rew v(ltom_incong_rew)];
incong_rew.htom = [htom_incong_rew v(htom_incong_rew)];

ltom_incong_opt(isnan(ltom_incong_opt)) = [];
htom_incong_opt(isnan(htom_incong_opt)) = [];

incong_opt.ltom = [ltom_incong_opt v(ltom_incong_opt)];
incong_opt.htom = [htom_incong_opt v(htom_incong_opt)];



ltom_cong_offer(isnan(ltom_cong_offer)) = [];
htom_cong_offer(isnan(htom_cong_offer)) = [];

cong_offer.ltom = [ltom_cong_offer v(ltom_cong_offer)];
cong_offer.htom = [htom_cong_offer v(htom_cong_offer)];

ltom_cong_rew(isnan(ltom_cong_rew)) = [];
htom_cong_rew(isnan(htom_cong_rew)) = [];

cong_rew.ltom = [ltom_cong_rew v(ltom_cong_rew)];
cong_rew.htom = [htom_cong_rew v(htom_cong_rew)];

ltom_cong_opt(isnan(ltom_cong_opt)) = [];
htom_cong_opt(isnan(htom_cong_opt)) = [];

cong_opt.ltom = [ltom_cong_opt v(ltom_cong_opt)];
cong_opt.htom = [htom_cong_opt v(htom_cong_opt)];

