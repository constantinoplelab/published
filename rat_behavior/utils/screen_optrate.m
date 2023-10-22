function [popt] = screen_optrate(S)
%calculate the probability of opting out for each session, ignoring violation trials.
  %inputs: S struct.
  %outputs: opt-out probability for each session.
  %modified by cmc 10/19/21.
  
popt = nan(length(S.pd),1);

for j = 1:length(S.pd)
    if isstruct(S.pd{j})
      goods = ~isnan(S.pd{j}.vios); %exclude violation trials.
      popt(j,1) = mean(S.pd{j}.optout(goods), 'omitnan')./mean(S.pd{j}.ProbCatch(goods), 'omitnan');
    end
end
