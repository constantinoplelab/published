function HJfig1d(datadir)
% Plot time to opt-out by block type for an example rat and all
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

Afolder = fullfile(datadir, 'data-published', 'A_structs');

rew = [5 10 20 40 80]; % list of rewards available (uL)
% initialize variables
high = nan(length(ratList), length(rew));
low = high;
mixed = high;

for r=1:length(ratList)
    % loop through ratlist to compute time to opt-out by block
    rat = upper(ratList{r});
    file = strcat('ratTrial_', rat);
    load(fullfile(Afolder, file), 'A');

    [hi, lo, mix] = blocks(A);

    high(r,:) = hi.wt(1,:);
    low(r,:) = lo.wt(1,:);
    mixed(r,:) = mix.wt(1,:);

    if strcmpi(ratList{r}, examprat)
        myer = [hi.er; lo.er; mix.er];
    end
end

x = 1:length(rew);

figure;
t = tiledlayout(1,2);
set(gcf, units='Inches', position=[6.5,4.4,6.75,3.73])

nexttile(1);
ex = find(cellfun(@(x)strcmpi(x, examprat), ratList));
shadedErrorBar(x, mixed(ex,:), myer(3,:), 'lineprops', ...
    {'color', [.41 .1 .81]}); hold on
shadedErrorBar(x, high(ex,:), myer(1,:), 'lineprops', '-r'); 
shadedErrorBar(x, low(ex,:), myer(2,:), 'lineprops', '-b');

xlim([0.5 5.5]);
ylim([9.5 15]);
yticks(10:15)
ylabel('Time to opt-out (s)');
set(gca, fontsize=15, tickdir='out', box='off', xtick=x, ...
    xticklabels={'5'; '10'; '20'; '40'; '80'});
subtitle('Example rat')

nexttile(2);
norm = mixed(:,3); % normalize by mixed 20uL
m = mixed ./ norm;
l = low ./ norm;
h = high ./ norm;
shadedErrorBar(x, mean(m, 'omitnan'), std(m, 'omitnan')./sqrt(length(ratList)), ...
    'lineprops', {'color', [.41 .1 .81]}); hold on
shadedErrorBar(x, mean(h, 'omitnan'), std(h, 'omitnan')./sqrt(length(ratList)), ...
    'lineprops', '-r');
shadedErrorBar(x, mean(l, 'omitnan'), std(l, 'omitnan')./sqrt(length(ratList)), ...
    'lineprops', '-b');

xlim([0.5 5.5]);
yticks(0.9:0.1:1.3)
ylabel('Fractional time to opt-out');
set(gca, fontsize=15, tickdir='out', box='off', xTick=x, ...
    XTickLabels={'5'; '10'; '20'; '40'; '80'});
subtitle(sprintf('N = %i rats', length(ratList)))



end