function plot_fig4f(datadir)
% Plot z-scored time to opt-out during mixed blocks conditioned on reward
% delay on the previous trial, averaged across rats. N = 16 rats.
% Short/long delay trials are defined the bottom/top quartile of the reward
% delay distribution across sessions for each rat.

% load all rats
load(fullfile(datadir, 'data-published/ratlist.mat'), 'ratList');

folder = fullfile(datadir, 'data-published/A_structs');
postshort = nan(length(ratList),5);
postlong = nan(length(ratList),5);
for rr=1:length(ratList)
    rat = ratList{rr};
    fprintf('%i of %i: %s\n', rr, length(ratList), rat);
    file = strcat(['ratTrial_', rat, '.mat']);
    load(fullfile(folder, file), 'A');
    [postshort(rr,:), postlong(rr,:)] = wt_vs_prevDelay(A);
end

% stat test
p = nan(1,5);
for col=1:5
    p(col) = signrank(postshort(:,col), postlong(:,col));
end
disp(p)

figure;
set(gcf, units='centimeters', renderer='painters', position=[24,10,6.95,7.64])
plotPretty(1:5, postshort, '#ECB8D1')
hold on
plotPretty(1:5, postlong, '#603A20')
set(gca, xtick=1:5, xticklabels={'5', '10', '20', '40', '80'}, ...
    ytick=-0.4:0.4:0.8, fontsize=12)
legend({'Post short', 'Post long'}, box='off', location='northwest')
ylabel('Time to opt-out (z-score)')
xlabel('Offered reward (uL)')
subtitle(sprintf('N = %i rats', length(ratList)))

end