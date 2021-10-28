function group = cvalgroup(Dat,type)
%designs labels for cross-validation groups for stratifying indices
%
%inputs:
%Dat: the data from OFC
%type: switch case for which kind of groups
%
%output:
%group: index for each sample that corresponds to a unique group


hits = Dat.hits;
prev_hits = [nan;hits(1:end-1)];

chosenprob = (Dat.went_right==1).*Dat.right_prob + ...
    (Dat.went_right==0).*Dat.left_prob;
chosenprob(isnan(Dat.went_right)) = nan;

chosenval = (Dat.went_right==1).*Dat.this_right_volume + ...
    (Dat.went_right==0).*Dat.this_right_volume;
chosenval(isnan(Dat.went_right)) = nan;
    
isRightSafe=  Dat.right_prob==1;
isLeftSafe=   Dat.left_prob==1;    
choseSafe = (isRightSafe & Dat.went_right==1) | ...
        (isLeftSafe & Dat.went_right==0);
    
XprevReward = prev_hits ==1  & ~isnan(hits);
Xreward = hits == 1;
Xsafe = choseSafe & ~isnan(hits); %consider change to hits==1.
Xplay =  ~isnan(hits);

ns = numel(Xplay);
group = cell(1,ns);

switch type
    case 'optout'
        group = Xplay;
    case 'behavior'
        %build group indices for combinations of
        %safe(1) vs risky(0)
        %rewared(1) vs loss(0)
        %prev reward(1) vs loss(0)
        %left (1) or right(0)
        %play or violate should be balanced themselves
        group_base = Xplay.*[Xsafe,Xreward,XprevReward,Dat.went_right==0];
        group_base = bi2de(group_base);
        group = group_base;
        
    case 'volume'
        %build labels for combinations of volume
        group_base = chosenval;
        group_base(isnan(group_base)) = 0;
        group = group_base;
        
    case 'prob'
        %build labels for different probabilities
        group_base = chosenprob;
        group_base(isnan(group)) = 0;
        for j = 1:ns
            group{j} = group_base(j);
        end
        
    case 'environment'
        %combo of 'volume' and 'prob'
        group_base = Xplay.*[chosenval,chosenprob];
        group = cell(1,ns);
        for j = 1:ns
            group{j} = num2str(group_base(j,:));
        end
end



