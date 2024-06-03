function HJfig1c(datadir)
% Plot time to initiate trial by block type for an example rat and all
% rats included in fiber photometry experiment (N = 17).
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

% load list of photometry rats and example rat
DArats = load(fullfile(datadir, 'data-published', 'ratlist_DA.mat'), ...
    'ratList');
AChrats = load(fullfile(datadir, 'data-published', 'ratlist_ACh.mat'), ...
    'ratList');
ratList = [DArats.ratList, AChrats.ratList];
examprat = 'J037';

Afolder = fullfile(datadir, 'data-published', 'A_Structs');

Lat = struct;
Lat.high.raw = nan(length(ratList), 1);
Lat.low.raw = nan(length(ratList), 1);
Lat.high.zscore = nan(length(ratList), 1);
Lat.low.zscore = nan(length(ratList), 1);

for rr=1:length(ratList)
    % loop through ratlist to compute trial initiation time by block
    rat = upper(ratList{rr});
    file = strcat('ratTrial_', rat);
    load(fullfile(Afolder, file), 'A');

    l = A.ITI;
    usethese = ~isnan(l);
    L = l(usethese);
    % filter outliers in trial initiation time
    outlier = [prctile(L,1), prctile(L,99)];
    L(L<outlier(1)) = NaN;
    L(L>outlier(2)) = NaN;
    blk = A.block(usethese);

    if strcmp(ratList{rr}, examprat)
        usethese_mi = blk==1;
        usethese_hi = blk==2;
        usethese_lo = blk==3;
        
        Lat.example_rat =...
            [mean(L(usethese_lo), 'omitnan'),...
            mean(L(usethese_mi), 'omitnan'),...
            mean(L(usethese_hi), 'omitnan');...
            std(L(usethese_lo), 'omitnan')./sqrt(sum(usethese_lo)),...
            std(L(usethese_mi), 'omitnan')./sqrt(sum(usethese_mi)), ...
            std(L(usethese_hi), 'omitnan')./sqrt(sum(usethese_hi))];
    end   
        
    Lat.high.raw(rr) = mean(L(blk==2), 'omitnan');
    Lat.low.raw(rr) = mean(L(blk==3), 'omitnan');
    Lat.mixed.raw(rr) = mean(L(blk==1), 'omitnan');
    
    L = (L-mean(L, 'omitnan'))./std(L, 'omitnan');
    Lat.high.zscore(rr) = mean(L(blk==2), 'omitnan');
    Lat.low.zscore(rr) = mean(L(blk==3), 'omitnan');
    Lat.mixed.zscore(rr) = mean(L(blk==1), 'omitnan');
end

figure
t = tiledlayout(1,2);
set(gcf, units='Inches', position=[5,4.7,7.76,3.73])

nexttile(1); hold on
errorbar(1, Lat.example_rat(1,1), Lat.example_rat(2,1),...
    'b_', markerfacecolor='b', linewidth=2, capsize=0)
errorbar(2, Lat.example_rat(1,2), Lat.example_rat(2,2),...
    '_', 'color', [.41 .1 .81], 'markerfacecolor', [.41 .1 .81], ...
    linewidth=2, capsize=0)
errorbar(3, Lat.example_rat(1,3), Lat.example_rat(2,3),...
    'r_', markerfacecolor='r', linewidth=2, capsize=0)

set(gca, fontsize=15, tickdir='out')
xlim([0.5 3.5])
ylim([2.6 6])
xticks(1:3)
xticklabels({'Low', 'Mixed', 'High'})
yticks(3:6)
ylabel('Time to initiate trial (s)')
subtitle('Example rat')

nexttile(2); hold on
errorbar(1, mean(Lat.low.zscore), std(Lat.low.zscore)/sqrt(length(ratList)),...
    'b_', markerfacecolor='b', linewidth=2, capsize=0)
errorbar(2, mean(Lat.mixed.zscore), std(Lat.mixed.zscore)/sqrt(length(ratList)),...
    '_', 'color', [.41 .1 .81], 'markerfacecolor', [.41 .1 .81], linewidth=2, capsize=0)
errorbar(3, mean(Lat.high.zscore), std(Lat.high.zscore)/sqrt(length(ratList)),...
    'r_', markerfacecolor='r', linewidth=2, capsize=0)

set(gca, fontsize=15, tickdir='out')
xlim([0.5 3.5])
xticks(1:3)
xticklabels({'Low', 'Mixed', 'High'})
yticks([-0.1, 0, 0.1])
ylabel({'Time to initiate trial', '(z-score)'})
subtitle(sprintf('N = %i rats', length(ratList)))


end
