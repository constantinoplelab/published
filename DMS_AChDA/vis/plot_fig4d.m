function plot_fig4d()

tau = 2.5; 
pcatches = [0.15, 0.2, 0.25];
dt = 0.1;
maxT = 30.; % approximately top 99 prctile given our task delay distribution
bins = 0.5:dt:maxT;
mycolor = {'#BCBEC0', '#6D6E71', '#231F20'};
kappa = 0.1;
D = 0.85;
R = [5,10,20,40,80];

figure;
set(gcf, units='inches', renderer='painters', ...
    position=[7,4,4.2,1.74])
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% First tile spans two columns (2x width)
ax1 = nexttile([1 2]);
hold on

% Second tile (1x width)
ax2 = nexttile;
hold on

for ix = 1:length(pcatches)
    p = pcatches(ix);
    C = 1 - p;
    V = log2(20)/tau * C * exp(-bins./tau) ./ (1 - C + C * exp(-bins./tau));
    WT = D * tau * log(C / (1 - C) * (R - kappa * tau) ./ (kappa * tau));

    plot(ax1, log2(bins), V, 'Color', mycolor{ix}, 'LineWidth', 1);
    plot(ax2, 1:5, WT, 'Color', mycolor{ix}, 'LineWidth', 1);
end

% Format first tile
set(ax1, 'TickDir', 'out', 'Box', 'off', 'FontSize', 8, ...
    'YTick', [0 1], 'XTick', 1:1:4);
xlim(ax1, [1 4]);
xlabel(ax1, 'Time in trial (s, log_{2})');
ylabel(ax1, 'V(t)');

% Format second tile
set(ax2, 'TickDir', 'out', 'Box', 'off', 'FontSize', 8, ...
    'XTick', 1:5, 'XTickLabel', {'5', '10', '20', '40', '80'});
xlim(ax2, [0.5 5.5]);
ylim(ax2, [7 17])
xlabel(ax2, 'Offered reward');
ylabel(ax2, 'Time to opt-out (s)');

end