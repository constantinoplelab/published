function plot_fig3c(datadir)
% Plot simulated value and trial initiation time from volume-modulated
% RPE model (left) and block-average trial initiation times overlaid with
% rat data (right).
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

% These model hyperparameters were determined by simulating with a grid of
% parameter combinations and selecting those that qualitatively best
% captured rats' data
alpha = 0.5; % learning rate
scalingfactor = 0.48; % scaling factor 

[modelITI, modelValue, ratTrial] = run_vanillaTD(datadir, 'reward', ...
    alpha, scalingfactor);

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
yl = [1 5];
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

% Plot average model simulated trial initiation time in low and high blocks
ratList = fields(ratTrial);
modelITI_low = nan(length(ratList),1);
modelITI_mixed = nan(length(ratList),1);
modelITI_high = nan(length(ratList),1);
ratITI_low = nan(length(ratList),1);
ratITI_mixed = nan(length(ratList),1);
ratITI_high = nan(length(ratList),1);

for r=1:length(ratList)
    ratname = ratList{r};
    Block = ratTrial.(ratname).block; % low block=3, high block=2
    modelITI_low(r) = mean(modelITI.(ratname)(Block==3), 'omitnan');
    modelITI_mixed(r) = mean(modelITI.(ratname)(Block==1), 'omitnan');
    modelITI_high(r) = mean(modelITI.(ratname)(Block==2), 'omitnan');
    ratITI_low(r) = mean(ratTrial.(ratname).ITI(Block==3), 'omitnan');
    ratITI_mixed(r) = mean(ratTrial.(ratname).ITI(Block==1), 'omitnan');
    ratITI_high(r) = mean(ratTrial.(ratname).ITI(Block==2), 'omitnan');
end

nexttile(inner, [1 1]);
errorbar([mean(modelITI_low), mean(modelITI_mixed), mean(modelITI_high)],...
    [sem(modelITI_low), sem(modelITI_mixed), sem(modelITI_high)], ...
    'k_', capsize=0); hold on
errorbar([mean(ratITI_low), mean(ratITI_mixed), mean(ratITI_high)], ...
    [sem(ratITI_low), sem(ratITI_mixed), sem(ratITI_high)], ...
    'b_', capsize=0)

xlim([.5 3.5]);
ylim([1.5 4])
set(gca, tickdir='out', box='off', xtick=[1 2 3], ...
    xticklabels={'Low', 'Mixed', 'High'}, ytick=2:4);
xlabel('Block type');
ylabel('Trial initiation time (s)')
legend({'model', 'rat'})

end

function [modelITI, modelValue, ratTrial, modelRPE] = ...
    run_vanillaTD(datadir, modeltype, alpha, scalingfactor)

% Load behavior data of all DA rats
ratList = loadRats(datadir, 'da');

modelITI = struct;
modelValue = struct;
ratTrial = struct;
modelRPE = struct;
for r=1:length(ratList)
    ratname = ratList{r};
    disp(ratname)
    fname = strcat('ratTrial_', ratname, '.mat');
    load(fullfile(datadir, 'data-published', 'A_structs', fname), 'A');
    [modelITI.(ratname), modelValue.(ratname), ratTrial.(ratname),...
        modelRPE.(ratname)] = ...
        vanillaTD_perRat(A, modeltype, alpha, scalingfactor);
end


end