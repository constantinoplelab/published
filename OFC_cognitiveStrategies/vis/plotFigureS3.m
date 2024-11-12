function plotFigureS3(incongruentPath_v2, incongruentPath_piriform, ...
    incongruentPath_motor, ephysPath)

%general
setDefaultFigProps

timeBin = [-4 8];
x = timeBin(1):0.05:timeBin(2);
wndw = [0 0.5]; %window for determining significance
w = arrayfun(@(y) find(x == wndw(y)), 1:2); %indices of window for computing significance

sem = @(x) std(x,'omitnan')./sqrt(sum(any(~isnan(x), 2)));

bonferroniCorrect = 2;

%% Block decoder

decoderPath_OFC = [ephysPath, filesep, 'Expert', filesep, 'Decoder', filesep];
decoderPath_naive = [ephysPath, filesep, 'Naive', filesep, 'Decoder', filesep];
decoderPath_v2 = [ephysPath, filesep, 'V2', filesep, 'Decoder', filesep];
decoderPath_pir = [ephysPath, filesep, 'Piriform', filesep, 'Decoder', filesep];
decoderPath_motor = [ephysPath, filesep, 'Motor', filesep, 'Decoder', filesep];

sem = @(x) std(x,'omitnan')./sqrt(sum(any(~isnan(x), 2)));

% block decoder performance for experts
matfiles = dir(fullfile(decoderPath_OFC, '*.mat'));
nsess = length(matfiles);

perfOFC = nan(nsess, 1);
for ii = 1:length(matfiles)
    load([matfiles(ii).folder filesep matfiles(ii).name], 'DecoderOutput');
    
    perfOFC(ii) = max(DecoderOutput.block.perf);
end

% block decoder performance for naive
matfiles = dir(fullfile(decoderPath_naive, '*.mat'));
nsess = length(matfiles);

perfNaive = nan(nsess, 1);
for ii = 1:length(matfiles)
    load([matfiles(ii).folder filesep matfiles(ii).name], 'DecoderOutput');
    
    perfNaive(ii) = max(DecoderOutput.block.perf);
end

% block decoder performance for v2

matfiles = dir(fullfile(decoderPath_v2, '*.mat'));
nsess = length(matfiles);

perfV2 = nan(nsess, 1);
for ii = 1:length(matfiles)
    load([matfiles(ii).folder filesep matfiles(ii).name], 'DecoderOutput');
    
    perfV2(ii) = max(DecoderOutput.block.perf);
end

% block decoder performance for piriform cortex

matfiles = dir(fullfile(decoderPath_pir, '*.mat'));
nsess = length(matfiles);

perfPir= nan(nsess, 1);
for ii = 1:length(matfiles)
    load([matfiles(ii).folder filesep matfiles(ii).name], 'DecoderOutput');
    
    perfPir(ii) = max(DecoderOutput.block.perf);
end

% block decoder performance for motor cortex

matfiles = dir(fullfile(decoderPath_motor, '*.mat'));
nsess = length(matfiles);

perfMotor= nan(nsess, 1);
for ii = 1:length(matfiles)
    load([matfiles(ii).folder filesep matfiles(ii).name], 'DecoderOutput');
    
    perfMotor(ii) = max(DecoderOutput.block.perf);
end


%% incongruent responses - v2
v2 = load([incongruentPath_v2, filesep, 'incongruentResponses.mat']);

incongPref_Avg_v2 = mean(v2.incongPref{4}, 1, 'omitnan');
incongPref_Sem_v2 = sem(v2.incongPref{4});
nIp_v2 = sum(any(~isnan(v2.incongPref{4}), 2));

incongNonpref_Avg_v2 = mean(v2.incongNonpref{4}, 1, 'omitnan');
incongNonpref_Sem_v2 = sem(v2.incongNonpref{4});
nInp_v2 = sum(any(~isnan(v2.incongNonpref{4}), 2));

congPref_Avg_v2 = mean(v2.congPref{4}, 1, 'omitnan');
congPref_Sem_v2 = sem(v2.congPref{4});
nCp_v2 = sum(any(~isnan(v2.congPref{4}), 2));

congNonpref_Avg_v2 = mean(v2.congNonpref{4}, 1, 'omitnan');
congNonpref_Sem_v2 = sem(v2.congNonpref{4});
nCnp_v2 = sum(any(~isnan(v2.congNonpref{4}), 2));


%stats
inNp_v2 = mean(v2.incongNonpref{4}(:, w(1):w(2)), 2, 'omitnan');
inNp_v2 = inNp_v2(~isnan(inNp_v2));

inP_v2 = mean(v2.incongPref{4}(:, w(1):w(2)), 2, 'omitnan');
inP_v2 = inP_v2(~isnan(inP_v2));
[incongP_v2, ~] = permute_test(inNp_v2, inP_v2, 1000); 

cNp_v2 = mean(v2.congNonpref{4}(:, w(1):w(2)), 2, 'omitnan');
cNp_v2 = cNp_v2(~isnan(cNp_v2));

cP_v2 = mean(v2.congPref{4}(:, w(1):w(2)), 2, 'omitnan');
cP_v2 = cP_v2(~isnan(cP_v2));
[congP_v2, ~] = permute_test(cNp_v2, cP_v2, 1000); 

%% incongruent responses - piriform
pir = load([incongruentPath_piriform, filesep, 'incongruentResponses.mat']);

incongPref_Avg_pir = mean(pir.incongPref{4}, 1, 'omitnan');
incongPref_Sem_pir = sem(pir.incongPref{4});
nIp_pir = sum(any(~isnan(pir.incongPref{4}), 2));

incongNonpref_Avg_pir = mean(pir.incongNonpref{4}, 1, 'omitnan');
incongNonpref_Sem_pir = sem(pir.incongNonpref{4});
nInp_pir = sum(any(~isnan(pir.incongNonpref{4}), 2));

congPref_Avg_pir = mean(pir.congPref{4}, 1, 'omitnan');
congPref_Sem_pir = sem(pir.congPref{4});
nCp_pir = sum(any(~isnan(pir.congPref{4}), 2));

congNonpref_Avg_pir = mean(pir.congNonpref{4}, 1, 'omitnan');
congNonpref_Sem_pir = sem(pir.congNonpref{4});
nCnp_pir = sum(any(~isnan(pir.congNonpref{4}), 2));


%stats
inNp_pir = mean(pir.incongPref{4}(:, w(1):w(2)), 2, 'omitnan');
inNp_pir = inNp_pir(~isnan(inNp_pir));

inP_pir = mean(pir.incongNonpref{4}(:, w(1):w(2)), 2, 'omitnan');
inP_pir = inP_pir(~isnan(inP_pir));
[incongP_pir, ~] = permute_test(inNp_pir, inP_pir, 1000); 

cNp_pir = mean(pir.congPref{4}(:, w(1):w(2)), 2, 'omitnan');
cNp_pir = cNp_pir(~isnan(cNp_pir));

cP_pir = mean(pir.congNonpref{4}(:, w(1):w(2)), 2, 'omitnan');
cP_pir = cP_pir(~isnan(cP_pir));
[congP_pir, ~] = permute_test(cNp_pir, cP_pir, 1000);

%% incongruent responses - motor cortex (any dorsal cells)
motor = load([incongruentPath_motor, filesep, 'incongruentResponses.mat']);

incongPref_Avg_motor = mean(motor.incongPref{4}, 1, 'omitnan');
incongPref_Sem_motor = sem(motor.incongPref{4});
nIp_motor = sum(any(~isnan(motor.incongPref{4}), 2));

incongNonpref_Avg_motor = mean(motor.incongNonpref{4}, 1, 'omitnan');
incongNonpref_Sem_motor = sem(motor.incongNonpref{4});
nInp_motor = sum(any(~isnan(motor.incongNonpref{4}), 2));

congPref_Avg_motor = mean(motor.congPref{4}, 1, 'omitnan');
congPref_Sem_motor = sem(motor.congPref{4});
nCp_motor = sum(any(~isnan(motor.congPref{4}), 2));

congNonpref_Avg_motor = mean(motor.congNonpref{4}, 1, 'omitnan');
congNonpref_Sem_motor = sem(motor.congNonpref{4});
nCnp_motor = sum(any(~isnan(motor.congNonpref{4}), 2));


%stats
inNp_motor = mean(motor.incongPref{4}(:, w(1):w(2)), 2, 'omitnan');
inNp_motor = inNp_motor(~isnan(inNp_motor));

inP_motor = mean(motor.incongNonpref{4}(:, w(1):w(2)), 2, 'omitnan');
inP_motor = inP_motor(~isnan(inP_motor));
[incongP_motor, ~] = permute_test(inNp_motor, inP_motor, 1000); 

cNp_motor = mean(motor.congPref{4}(:, w(1):w(2)), 2, 'omitnan');
cNp_motor = cNp_motor(~isnan(cNp_motor));

cP_motor = mean(motor.congNonpref{4}(:, w(1):w(2)), 2, 'omitnan');
cP_motor = cP_motor(~isnan(cP_motor));
[congP_motor, ~] = permute_test(cNp_motor, cP_motor, 1000); 

%% 
% plot 
all = {perfOFC perfNaive perfV2 perfPir perfMotor};
xvec = 1:length(all);

jitx = arrayfun(@(x) randn(length(all{x}), 1) ./ 25 + xvec(x), xvec, ...
    'UniformOutput', false);

groups = cell2mat(arrayfun(@(x) ones(length(all{x}),1)*xvec(x), xvec, ...
    'UniformOutput', false)');


figure; hold on
t = tiledlayout(4, 2, 'TileSpacing', 'compact');

nexttile([1,2])
hold on
arrayfun(@(x) scatter(jitx{x}, all{x}, 20, 'k', 'filled', ...
    'markerfacealpha', 0.15), xvec)
hold on
plot([0.5 5.5], [0.33 0.33], '--k')
arrayfun(@(x) errorbar(xvec(x)+0.25, mean(all{x}), std(all{x}), '_k'), xvec)
xlim([0.5 5.5])
ylim([0 1])
xticks(xvec)
yticks(0:0.2:1)
ylim([0 1])
xticklabels({'Expert' 'Naive' 'V2' 'Pir' 'Motor'})
ylabel('Decoder performance')
set(gca, 'TickDir', 'out'); box off;

ax1 = nexttile;
shadedErrorBar(x, incongNonpref_Avg_v2, incongNonpref_Sem_v2, 'lineprops', ...
    {'color', 'k', 'linewidth', 1}) %preferred transition
shadedErrorBar(x, incongPref_Avg_v2, incongPref_Sem_v2, 'lineprops', ...
    {'color', [0.65 0.65 0.65], 'linewidth', 1}) %non-preferred transition
hold on
plot([0 0], [0 25], '--k')
text(2, 1, string(incongP_v2*bonferroniCorrect))
yticks(0:5:25)
set(gca, 'TickDir', 'out'); box off;
axis square
ax1.YRuler.TickLabelGapOffset = 1;
title('\rm Incongruent')

ax2 = nexttile;
shadedErrorBar(x, congNonpref_Avg_v2, congNonpref_Sem_v2, 'lineprops', ...
    {'color', 'k', 'linewidth', 1})
shadedErrorBar(x, congPref_Avg_v2, congPref_Sem_v2, 'lineprops', ...
    {'color', [0.65 0.65 0.65], 'linewidth', 1})
hold on
plot([0 0], [0 25], '--k')
text(2, 1, string(congP_v2*bonferroniCorrect))
yticks(0:5:25)
set(gca, 'TickDir', 'out'); box off;
axis square
ax2.YRuler.TickLabelGapOffset = 1;
title('\rm Congruent') 
linkaxes([ax1 ax2])
xlim([-1 3])
ylim([0 25])

ax3 = nexttile;
shadedErrorBar(x, incongNonpref_Avg_pir, incongNonpref_Sem_pir, 'lineprops', ...
    {'color', 'k', 'linewidth', 1}) %preferred transition
shadedErrorBar(x, incongPref_Avg_pir, incongPref_Sem_pir, 'lineprops', ...
    {'color', [0.65 0.65 0.65], 'linewidth', 1}) %non-preferred transition
hold on
plot([0 0], [0 15], '--k')
text(2, 1, string(incongP_pir*bonferroniCorrect))
yticks(0:5:15)
set(gca, 'TickDir', 'out'); box off;
axis square
ax3.YRuler.TickLabelGapOffset = 1;
legend('Pref', 'Nonpref', 'Location', 'northeastoutside')

ax4 = nexttile;
shadedErrorBar(x, congNonpref_Avg_pir, congNonpref_Sem_pir, 'lineprops', ...
    {'color', 'k', 'linewidth', 1})
shadedErrorBar(x, congPref_Avg_pir, congPref_Sem_pir, 'lineprops', ...
    {'color', [0.65 0.65 0.65], 'linewidth', 1})
hold on
plot([0 0], [0 15], '--k')
text(2, 1, string(congP_pir*bonferroniCorrect))
yticks(0:5:15)
set(gca, 'TickDir', 'out'); box off;
axis square
ax4.YRuler.TickLabelGapOffset = 1;
linkaxes([ax3 ax4])
xlim([-1 3])
ylim([0 15])

ax5 = nexttile;
shadedErrorBar(x, incongNonpref_Avg_motor, incongNonpref_Sem_motor, 'lineprops', ...
    {'color', 'k', 'linewidth', 1}) %preferred transition
shadedErrorBar(x, incongPref_Avg_motor, incongPref_Sem_motor, 'lineprops', ...
    {'color', [0.65 0.65 0.65], 'linewidth', 1}) %non-preferred transition
hold on
plot([0 0], [0 15], '--k')
text(2, 1, string(incongP_motor*bonferroniCorrect))
yticks(0:5:15)
set(gca, 'TickDir', 'out'); box off;
axis square
ax5.YRuler.TickLabelGapOffset = 1;

ax6 = nexttile;
shadedErrorBar(x, congNonpref_Avg_motor, congNonpref_Sem_motor, 'lineprops', ...
    {'color', 'k', 'linewidth', 1})
shadedErrorBar(x, congPref_Avg_motor, congPref_Sem_motor, 'lineprops', ...
    {'color', [0.65 0.65 0.65], 'linewidth', 1})
hold on
plot([0 0], [0 15], '--k')
text(2, 1, string(congP_motor*bonferroniCorrect))
yticks(0:5:15)
set(gca, 'TickDir', 'out'); box off;
axis square
ax6.YRuler.TickLabelGapOffset = 1;
linkaxes([ax5 ax6])
xlim([-1 3])
ylim([0 15])
ylabel(ax3, 'Firing rate (Hz)', 'FontSize', 12)
xlabel(t, 'Time from reward (s)')
set(gcf,'renderer','painter')
set(gcf, 'units', 'centimeters', 'position', [10 5 11 17])
