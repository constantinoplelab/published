function plot_fig3f(datadir)
% Plot simulated value and trial initiation time from delay-modulated
% RPE model (left) and block-average trial initiation times overlaid with
% rat data (right).
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

% These model hyperparameters were determined by simulating with a grid of
% parameter combinations and selecting those that qualitatively best
% captured rats' data
alpha = 0.46; % learning rate
scalingfactor = 0.54; % scaling factor

[modelITI, modelValue, ratTrial] = run_vanillaTD(datadir, 'delay', alpha, ...
    scalingfactor);

exampleRat = 'J029';
blkLength = 40; % each block is 40 trials

% Plot
figure;
set(gcf, units='inches', position=[6,5,7,2.5])
main = tiledlayout(1,1);

inner = tiledlayout(main, 1,3, tilespacing = 'compact', ...
    padding='compact');

% Example rat
% Plot model simulated value and trial initiation time for 3 consecutive
% blocks
nexttile(inner, [1 2]);
yl = [0 12];
mixedcolor=[0.4823529411764706 0.1568627450980392 0.48627450980392156];
fill([0 blkLength blkLength 0], [yl(1) yl(1) yl(2) yl(2)],...
    mixedcolor, facealpha=0.15, edgecolor='none'); hold on
fill([blkLength blkLength*2 blkLength*2 blkLength], ...
    [yl(1) yl(1) yl(2) yl(2)], 'b', facealpha=0.15, edgecolor='none');
fill([blkLength*2 blkLength*3 blkLength*3 blkLength*2], ...
    [yl(1) yl(1) yl(2) yl(2)], mixedcolor, facealpha=0.15, edgecolor='none');
fill([blkLength*3 blkLength*4 blkLength*4 blkLength*3], ...
    [yl(1) yl(1) yl(2) yl(2)], 'r', facealpha=0.15, edgecolor='none');


V_ex = modelValue.(exampleRat);
ITI_ex = modelITI.(exampleRat);
% Exclude violation trials
V_ex = V_ex(ratTrial.(exampleRat).vios~=1);
ITI_ex = ITI_ex(ratTrial.(exampleRat).vios~=1);
% fill in NaNs with a moving average filter
V_ex = fillmissing(V_ex, 'movmean', 5);
ITI_ex = fillmissing(ITI_ex, 'movmean', 5);

yyaxis left
plot(V_ex(1:blkLength*4), '-', color='k', linewidth=0.5); 
hold on
xlim([blkLength blkLength*4])
for xx=blkLength:blkLength:blkLength*4
    plot([xx, xx+blkLength], ...
        [mean(V_ex(xx:xx+blkLength)), mean(V_ex(xx:xx+blkLength))], ...
        'k-', linewidth=1.5);
end
ylim(yl)
yticks([])
xticks([])
ylabel('Value')

yyaxis right
plot(ITI_ex(1:blkLength*4), '-', color=[0 0 0 0.25], linewidth=0.5);
hold on
for xx=blkLength:blkLength:blkLength*4
    plot([xx, xx+blkLength], ...
        [mean(ITI_ex(xx:xx+blkLength)), mean(ITI_ex(xx:xx+blkLength))], ...
        '-', color=[0 0 0 0.25], linewidth=1.5);
end
yticks([])
xlabel('Trial')
ylabel('Time to initiate trial')
set(gca, tickdir='out', box='off');

% Plot average model simulated trial initiation time conditioned on
% previous delay
ratList = fields(ratTrial);
modelITI_short = nan(length(ratList),1);
modelITI_long = nan(length(ratList),1);
ratITI_short = nan(length(ratList),1);
ratITI_long = nan(length(ratList),1);

for r=1:length(ratList)
    ratname = ratList{r};
    
    prevDelay = [nan; ratTrial.(ratname).reward_delay(1:end-1)];
    prevDelay(prevDelay==100) = nan;
    posthit = [0; ratTrial.(ratname).hits(1:end-1)];
    q = quantile(prevDelay, [0.25 0.5 0.75]);
    delay_bin = discretize(prevDelay, [-inf q inf]);
    
    usethese = posthit;
    modelITI_short(r) = mean(modelITI.(ratname)(delay_bin==1 & usethese), 'omitnan');
    modelITI_long(r) = mean(modelITI.(ratname)(delay_bin==4 & usethese), 'omitnan');
    ratITI_short(r) = mean(ratTrial.(ratname).ITI(delay_bin==1 & usethese), 'omitnan');
    ratITI_long(r) = mean(ratTrial.(ratname).ITI(delay_bin==4 & usethese), 'omitnan');
end

nexttile(inner, [1 1]);
errorbar([mean(modelITI_short), mean(modelITI_long)],...
    [sem(modelITI_short), sem(modelITI_long)], ...
    'k_', capsize=0); hold on
errorbar([mean(ratITI_short), mean(ratITI_long)], ...
    [sem(ratITI_short), sem(ratITI_long)], ...
    'b_', capsize=0)

xlim([.5 2.5]);
ylim([1.2 3.7])
set(gca, tickdir='out', box='off', xtick=[1 2], ...
    xticklabels={'Short', 'Long'}, ytick=2:3);
xlabel('Previous delay');
ylabel('Trial initiation time (s)')
legend({'model', 'rat'})

end