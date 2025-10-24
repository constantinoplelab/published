function plot_fig5ef(datadir)
% Plot fraction of cells with significant slopes at different task events
% (left) and histogram of slopes for significantly modulated cells at the
% offer cue (right).
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

fpath = fullfile(datadir, 'data-published/NpxData/pD1pD2.mat');
load(fpath, 'pval', 'slope')

hasData = ~isnan(pval(1,:));

figure
t = tiledlayout(1,2, padding='compact');
set(gcf, 'units', 'inches', 'position', [7,4,4.7,2])
pci = nan(sum(hasData),2); 

nexttile(2)
usethese = hasData & pval(1,:)<0.05 & slope(1,:)~=0;
H = histogram(slope(1,usethese), facecolor='none', binwidth=0.2);
hEdges = H.BinEdges;
histogram(slope(1, usethese&slope(1,:)>0), BinEdges=hEdges, ...
    facecolor='#004489')
hold on
histogram(slope(1, usethese&slope(1,:)<0), BinEdges=hEdges, ...
    facecolor='#b92237')
xline(0, 'k', linewidth=1.5)
ylabel('N (cell)')
xlabel('slope')
xlim([-4 4])
ylim([0 20])
set(gca, tickdir='out', box='off')

nexttile(1)
for a=1:size(slope,1)
    nsig = length(slope(a,pval(a,:)<0.05));
    frac = 100*nsig/sum(hasData);
    [~,pci(a,:)] = binofit(length(slope(a,pval(a,:)<0.05)), sum(hasData));
    bar(a, frac, facecolor='none', barwidth=0.35); hold on
    errorbar(a, frac, pci(a,1), pci(a,2), 'k')
end
xlim([0.5 3.5])
ylim([0 60])
set(gca, tickdir='out', box='off', xtick=1:3, ...
    xticklabels={'OfferCue', 'RewPortOn', 'RewCue'})
ylabel('Fraction significant (%)')

sgtitle(t, sprintf('N = %d cells', sum(hasData)), fontsize=9)


end