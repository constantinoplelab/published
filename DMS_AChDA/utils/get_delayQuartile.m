function quartiles = get_delayQuartile(datadir, ratname)

% get the bottom and top quartiles of reward delay to set lower and upper
% bounds
load(fullfile(datadir, 'data-published\A_structs', ...
    strcat('ratTrial_', ratname, '.mat')), 'A')

delays = A.reward_delay;
delays(delays==100) = [];
delays(delays<0.75) = [];
quartiles = quantile(delays, [0:1/4:1]);

end