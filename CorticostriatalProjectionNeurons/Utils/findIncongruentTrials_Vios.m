function [t_incongruent, t_congruent] = findIncongruentTrials_Vio(S)
%Finds the index for the first REWARDED incongruent trial in a mix block
%and the most recent rewarded trial congruent with the adaptation block.

%input: 
    % S: S struct of neuropixels data

%output: 
    % t_incongruent:
        % ltom: low to mix incongruent trial [trial, volume rewarded]
        % htom: high to mix
    %t_congruent: 
        % ltom: first congruent trial after the true transition into the 
        %       mixed block but are before the incongruent trial
        %       [trial, rewarded volume]
        % htom: high to mix
%reward volume = 1:5, to account for male and female rats

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

for tt = 1:length(htom) %for each high to mix transition
    next = [allChange(allChange > htom(tt)); length(v)];    % next = [block changes that do h2m; number of trials where there are rewards]
    
    i = intersect(lowRew, htom(tt):next(1));    %find trials where there are low volume rewards AND trials after the htomix block changes
    rewarded = i;       % find trials index of low vol trials, (hit, vio, catch)
    %rewarded = intersect(i, find(S.hits));  %find the trial numbers of low rewarded trials (1st should be the 1st rewarded incongruent trial)
    
    try %find trials congruent with the last adaption block
        c = intersect([highRew; twenty], 1:rewarded(1));     
    catch
        c = htom(tt):next(1); %no congruent trials
    end
    %c_rew = intersect(c, find(S.hits));
        c_rew = c;


    if ~isempty(rewarded)
        htom_incong_rew(tt) = rewarded(1);  % keep only the first incongruent trial (rewarded)
    else
        htom_incong_rew(tt) = nan;
    end

    if ~isempty(c_rew)
        htom_cong_rew(tt) = c_rew(end); %keep only congruent trial (rewarded) most recent to the incongruent trial
    else
        htom_cong_rew(tt) = nan;
    end
end

for tt = 1:length(ltom)
    next = [allChange(allChange > ltom(tt)); length(v)];

    i = intersect(highRew, ltom(tt):next(1));
    rewarded = i;   
    %rewarded = intersect(i, find(S.hits));
    opt = intersect(i, find(S.optout));
    
    try
        c = intersect([lowRew; twenty], 1:rewarded(1));
    catch
        c = ltom(tt):next(1); %if no incongruent trials
    end
    %c_rew = intersect(c, find(S.hits));
        c_rew = c;


    if ~isempty(rewarded)
        ltom_incong_rew(tt) = rewarded(1);
    else
        ltom_incong_rew(tt) = nan;
    end

    if ~isempty(c_rew)
        ltom_cong_rew(tt) = c_rew(end);
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
