function plotFigureS6_otherAreas(ephysPath, v2Path, motorPath, piriformPath, codePath)
% Plot figure S6. YOU MUST RUN processEphysData_experts FOR V2, MOTOR, AND
% PIRIFORM BEFORE RUNNING THIS FUNCITON. Save output of each run to a
% separate folder.

% Plot decoder performance for OFC, V2, Piriform cortex,
% and Motor cortex in experts. Plot responses for 20ul trials pre- and
% post- incongruent for non-OFC recording areas

% INPUTS:
%   ephysPath = local path to ephys data with decoder folders for each
%       recording area
%   v2Path = local path to the saved output from running
%       processEphysData_experts for V2 data
%   motorPath = local path to the saved output from running
%       processEphysData_experts for motor cortex data
%   piriformPath = local path to the saved output from running
%       processEphysData_experts for piriform cortex data
%   codePath = local path to code downloaded from github



%% Set up paths
s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codePath, 'IgnoreCase', ispc);

if ~onPath % only add code dir to path if it already isn't
    addpath(genpath(codePath))
end

%% general
timeBin = [-4 8]; % time window around each event saved in the SU structs
x = timeBin(1):0.05:timeBin(2);
wndw = [-0.5 1];
w = arrayfun(@(y) find(x == wndw(y)), 1:2); %indices of window for computing significance
xvec = w(1):w(2);

sem = @(x) std(x,'omitnan')./sqrt(sum(any(~isnan(x), 2)));

setDefaultFigProps

%% Block decoder

decoderPath_OFC = [ephysPath, filesep, 'Expert', filesep, 'OFC', filesep, ...
    'Decoder', filesep];
decoderPath_v2 = [ephysPath, filesep, 'Expert', filesep, 'V2', filesep, ...
    'Decoder', filesep];
decoderPath_pir = [ephysPath, filesep, 'Expert', filesep, 'Piriform', ...
    filesep, 'Decoder', filesep];
decoderPath_motor = [ephysPath, filesep, 'Expert', filesep, 'Motor', ...
    filesep, 'Decoder', filesep];

% block decoder performance for OFC
matfiles = dir(fullfile(decoderPath_OFC, '*.mat'));
nsess = length(matfiles);

perfOFC = nan(nsess, 1);
for jj = 1:length(matfiles)
    load([matfiles(jj).folder filesep matfiles(jj).name], 'DecoderOutput');
    
    perfOFC(jj) = max(DecoderOutput.block.perf);
end

% block decoder performance for v2

matfiles = dir(fullfile(decoderPath_v2, '*.mat'));
nsess = length(matfiles);

perfV2 = nan(nsess, 1);
for jj = 1:length(matfiles)
    load([matfiles(jj).folder filesep matfiles(jj).name], 'DecoderOutput');
    
    perfV2(jj) = max(DecoderOutput.block.perf);
end

% block decoder performance for piriform cortex

matfiles = dir(fullfile(decoderPath_pir, '*.mat'));
nsess = length(matfiles);

perfPir= nan(nsess, 1);
for jj = 1:length(matfiles)
    load([matfiles(jj).folder filesep matfiles(jj).name], 'DecoderOutput');
    
    perfPir(jj) = max(DecoderOutput.block.perf);
end

% block decoder performance for motor cortex

matfiles = dir(fullfile(decoderPath_motor, '*.mat'));
nsess = length(matfiles);

perfMotor= nan(nsess, 1);
for jj = 1:length(matfiles)
    load([matfiles(jj).folder filesep matfiles(jj).name], 'DecoderOutput');
    
    perfMotor(jj) = max(DecoderOutput.block.perf);
end


%% plot decoder panel
all = {perfOFC perfV2 perfPir perfMotor};
x2 = 1:length(all);

jitx = arrayfun(@(x) randn(length(all{x}), 1) ./ 25 + x2(x), x2, ...
    'UniformOutput', false);

% groups = cell2mat(arrayfun(@(x) ones(length(all{x}),1)*x2(x), x2, ...
%     'UniformOutput', false)');


figure; hold on
arrayfun(@(x) scatter(jitx{x}, all{x}, 20, 'k', 'filled', ...
    'markerfacealpha', 0.15), x2)
hold on
plot([0.5 5.5], [0.33 0.33], '--k')
arrayfun(@(x) errorbar(x2(x)+0.25, mean(all{x}), std(all{x}), '_k'), x2)
xlim([0.5 5.5])
ylim([0 1])
xticks(xvec)
yticks(0:0.2:1)
ylim([0 1])
xticklabels({'OFC' 'V2' 'Pir' 'Motor'})
ylabel('Decoder performance')
set(gca, 'TickDir', 'out'); box off;

set(gcf, 'units', 'centimeters', 'position', [10 10 10 4])
set(gcf,'renderer','painter')


%% plot incongruent responses for each area

paths = {v2Path motorPath piriformPath};
allYlims = {[-1 1.3] [-0.3 0.82] [-0.3 0.5]};

for ii = 1:length(paths)
    % load incongruent responses
    incongruentResponses = load([paths{ii}, 'incongruentResponses.mat']);

    % averages
    nAlign = 4;

    prePref_Avg = arrayfun(@(x) mean(incongruentResponses.pre_pref{x}, 1, 'omitnan'), ...
        1:nAlign, 'UniformOutput', false);
    prePref_Sem = arrayfun(@(x) sem(incongruentResponses.pre_pref{x}), 1:nAlign, ...
        'uniformoutput', false);

    preNonpref_Avg = arrayfun(@(x) mean(incongruentResponses.pre_nonpref{x}, 1, 'omitnan'), ...
        1:nAlign, 'UniformOutput', false);
    preNonpref_Sem = arrayfun(@(x) sem(incongruentResponses.pre_nonpref{x}), 1:nAlign, ...
        'uniformoutput', false);

    postPref_Avg = arrayfun(@(x) mean(incongruentResponses.post_pref{x}, 1, 'omitnan'), ...
        1:nAlign, 'UniformOutput', false);
    postPref_Sem = arrayfun(@(x) sem(incongruentResponses.post_pref{x}), 1:nAlign, ...
        'uniformoutput', false);

    postNonpref_Avg = arrayfun(@(x) mean(incongruentResponses.post_nonpref{x}, 1, 'omitnan'), ...
        1:nAlign, 'UniformOutput', false);
    postNonpref_Sem = arrayfun(@(x) sem(incongruentResponses.post_nonpref{x}), 1:nAlign, ...
        'uniformoutput', false);

    % compute p-values for each time bin
    pExp_post = cell(1, nAlign);
    pExp_pre = pExp_post;
    for jj = 1:nAlign
        pExp_post{jj} = arrayfun(@(y) permute_test(...
            incongruentResponses.post_nonpref{jj}(:,y), ...
            incongruentResponses.post_pref{jj}(:,y), 1000, 1), xvec);
        pExp_pre{jj} = arrayfun(@(y) permute_test(...
            incongruentResponses.pre_nonpref{jj}(:,y), ...
            incongruentResponses.pre_pref{jj}(:,y), 1000, 1), xvec);
    end

    % find significant bins
    pExp_pre_sig = arrayfun(@(x) pExp_pre{x} < 0.05, 1:nAlign, ...
        'uniformoutput', false);
    pExp_post_sig = arrayfun(@(x) pExp_post{x} < 0.05, 1:nAlign, ...
        'uniformoutput', false);
    

    figure; hold on
    tiledlayout(2, nAlign, 'TileSpacing', 'compact')

    %pre-incongruent
    labels = {'Trial start' 'LED on' 'LED off' 'Reward'};
    for jj = 1:nAlign

        sig = false(1, length(xvec));
        %find 3 consecutive significant bins
        a0 = [0 pExp_pre_sig{jj} 0]; % add 0 at beginning and end to find pattern. also shifts index so start is the first true value
        startTrue = strfind(a0,[0 1]);
        endTrue = strfind(a0,[1 0]);
        numSigBins = endTrue - startTrue;
        use = find(numSigBins > 2);
        good = arrayfun(@(x) startTrue(use(x)):endTrue(use(x))-1, 1:length(use), ...
            'uniformoutput', false);
        sig(cell2mat(good)) = true;

        % plot
        ax(jj) = nexttile;
        shadedErrorBar(x(xvec), prePref_Avg{jj}(xvec), prePref_Sem{jj}(xvec), ...
            'lineprops', {'color', [0.65 0.65 0.65], 'linewidth', 1})
        hold on
        shadedErrorBar(x(xvec), preNonpref_Avg{jj}(xvec), preNonpref_Sem{jj}(xvec), ...
            'lineprops', {'color', 'k', 'linewidth', 1})
        scatter(x(xvec(sig)), allYlims{ii}(2), '*k')
        xline(0, '--', 'linewidth', 0.5, 'color', [0.7 0.7 0.7])
        set(gca, 'TickDir', 'out'); box off;
        axis square
        ax(jj).YRuler.TickLabelGapOffset = 1;
        title(['\rm' labels{jj}])
    end
    ylabel(ax(1), 'Firing rate (Hz)')
    linkaxes(ax)
    xlim([-0.5 1])
    ylim(allYlims{ii})


    for jj = 1:nAlign

        sig = false(1, length(xvec));
        %find 3 consecutive significant bins
        a0 = [0 pExp_post_sig{jj} 0]; % add 0 at beginning and end to find pattern. also shifts index so start is the first true value
        startTrue = strfind(a0,[0 1]);
        endTrue = strfind(a0,[1 0]);
        numSigBins = endTrue - startTrue;
        use = find(numSigBins > 2);
        good = arrayfun(@(x) startTrue(use(x)):endTrue(use(x))-1, 1:length(use), ...
            'uniformoutput', false);
        sig(cell2mat(good)) = true;

        ax(jj) = nexttile;
        shadedErrorBar(x(xvec), postPref_Avg{jj}(xvec), postPref_Sem{jj}(xvec), ...
            'lineprops', {'color', [0.65 0.65 0.65], 'linewidth', 1})
        hold on
        shadedErrorBar(x(xvec), postNonpref_Avg{jj}(xvec), postNonpref_Sem{jj}(xvec), ...
            'lineprops', {'color', 'k', 'linewidth', 1})
        scatter(x(xvec(sig)), allYlims{ii}(2), '*k')
        xline(0, '--', 'linewidth', 0.5, 'color', [0.7 0.7 0.7])
        set(gca, 'TickDir', 'out'); box off;
        axis square
        ax(jj).YRuler.TickLabelGapOffset = 1;
        title(['\rm' labels{jj}])
    end
    ylabel(ax(1), 'Firing rate (Hz)')
    linkaxes(ax)
    xlim([-0.5 1])
    ylim(allYlims{ii})

    set(gcf, 'units', 'centimeters', 'position', [10 10 12 8])
    set(gcf,'renderer','painter')

end

