function [postshort, postlong, pval] = wt_vs_prevDelay(A, prevrew)
% Get wait time for mixed block 20uL by whether the previous trial reward
% was delivered after a long vs short delay (bottom vs top tercile)
% rew: conditioning on previous reward volume (1,2,3,4,5); default is all

%convert rewards to ordinal (log) scale
[A.reward, ~] = convertreward(A.reward);

postshort = nan(1,5); 
postlong = postshort; 
pval = postshort;

% get the bottom and top quartile of reward delay to set lower and upper
% bounds
delays = A.reward_delay;
delays(delays==100) = nan;
delays(delays<0.75) = nan;
quartiles = quantile(delays, [0.25 0.5 0.75]);
bins = discretize(delays, [-inf, quartiles, inf]);

posthit = [0; A.hits(1:end-1)];
prevBin = [nan; bins(1:end-1)];
prevR = [nan; A.reward(1:end-1)];

wt = A.wait_time;
wt(A.optout==0) = nan; % only get wait time on opt-out trials
wt = (wt-mean(wt, 'omitnan'))./std(wt, 'omitnan'); %z-score

for rew=1:5
    if nargin<2
        t1 = find(A.block==1 & A.reward==rew & prevBin==1 & A.optout==1 & ...
            posthit);
        t2 = find(A.block==1 & A.reward==rew & prevBin==4 & A.optout==1 & ...
            posthit);
    else
        t1 = find(A.block==1 & A.reward==rew & prevBin==1 & A.optout==1 & ...
            posthit & prevR==prevrew);
        t2 = find(A.block==1 & A.reward==rew & prevBin==4 & A.optout==1 & ...
            posthit & prevR==prevrew);
    end
    if ~isnan(mean(wt(t1), 'omitnan')) && ~isnan(mean(wt(t2), 'omitnan'))
        postshort(rew) = mean(wt(t1), 'omitnan');
        postlong(rew) = mean(wt(t2), 'omitnan');
        [~,pval(rew)] = kstest2(wt(t2), wt(t1));
    end
end

end