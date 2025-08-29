function plotFigureS3_muscimolPhys(muscimolPhysPath)
% Plot example cell recorded before and after muscimol infusion in an
% anesthetized rat. Plot summary of activity for all cells recorded before
% and after muscimol infusion.

% INPUTS:
%   muscimolPhysPath = local path to muscimol physiology data downloaded
%       from zenodo

%% process muscimol physiology data
out = muscimolPhysiology(muscimolPhysPath);

%% Plot  

figure; 
tiledlayout(1, 2, 'TileSpacing', 'compact')

%example repressed cell during muscimol infusion
idx = find(out.xvec == -15):find(out.xvec == 30);

nexttile
plot(out.xvec(idx), out.hmat(idx), 'k', 'linewidth', 1.5)
hold on
plot([0 0], [0 10], '--k', 'linewidth', 1)
xlabel({'Time from infusion'; '(minutes)'})
ylabel('Firing rate (Hz)')
xlim([-15 30])
set(gca, 'TickDir', 'out'); box off;
ax1 = gca;
ax1.YRuler.TickLabelGapOffset = 1;
title('\rm Example cell')
axis square

%muscimol physiology summary
nexttile
plot(out.x, out.sigFit, 'color', 'm', 'linewidth', 1.5)
hold on
scatter(out.bins, out.binnedratio, 100, '.k')
errorbar(out.bins, out.binnedratio, out.binnedEr, 'k', 'linestyle', 'none')
xl = get(gca, 'xlim');
plot([0 xl(2)], [1 1], '--k')
xlabel({'Distance from infusion'; 'site (mm)'})
ylabel('Normalized activity')
xlim([0 3.5])
ylim([0 2.5])
xticks(0:0.5:3.5)
text(0.2, 2.2, 'Sigmoid fit', 'color', 'm', 'FontSize', 8)
set(gca, 'TickDir', 'out'); box off;
ax2 = gca;
ax2.YRuler.TickLabelGapOffset = 1;
axis square

set(gcf,'renderer','painter')
set(gcf, 'units', 'centimeters', 'position', [10 10 8 6])


