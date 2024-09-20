function ProcessData_Figure2(BestFitEarly, BestFitLate, codedir, savedir)
%ProcessData_Figure2 - function to process data to generate Figure 2. Saves
% data in savedir to be accessed by PlotFigure2
% INPUTS:
% BestFitEarly - BestFit.mat file found at: Mah_CellReports_BehavorialData/FitAll_ITI_VanillaAlpha_First10/BestFit.mat
%   from Zenodo
% BestFitLate - BestFit.mat file found at: Mah_CellReports_BehavorialData/FitAll_ITI_VanillaAlpha_Last10/BestFit.mat
%   from Zenodo
% codedir - Directory of code (e.g., published/dynamic_learning_rate/)
% savedir - Location you would like the outputs to be saved

%% Set up paths

s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codedir, 'IgnoreCase', ispc);
if ~onPath % only add code dir to path if it already isn't
    addpath(genpath(codedir))
end

%% Process data

% Pull list of rats 
fieldlist = fields(BestFitEarly);
ratList = fieldlist(structfun(@isstruct, BestFitEarly));

N = 10; % Cutoff for early and late trials (This pulls first and last 10 trials) 
nback = 9; % Number of previous trials to regress over
twin = 40; % Window around block transition to pull for dynamics plot
smoothfactor = 10; % Size of smoothing window for dynamics plot
exampleRatMdl = 'J076'; % Example rat ID

% Preallocate data structures
[ltom, htom] = deal(nan(length(ratList), 2*twin+1));

[betasEarly, betasLate,...
    betasEarlySim, betasLateSim] = deal(nan(length(ratList), nback+1));
[tauEarly, tauLate] = deal(nan(size(ratList)));

[nLLEarlyEarly, nLLEarlyLate, nLLLateEarly, nLLLateLate] =...
    deal(nan(size(ratList)));

[expParamsEarly, expParamsLate] = deal(nan(length(ratList), 2));

% Loop over rats
for rr = 1:length(ratList)
    fprintf('%d out of %d\n', rr, length(ratList)) % Display progress

    % Pull data for rat (A) and for model (ALate)
    [A, ALate] = deal(BestFitEarly.(ratList{rr}).All.ratTrial);
    ALate.ITI = BestFitLate.(ratList{rr}).All.LatMdl; % Pull model ITI

    % Logical vectors to pull early and late trials
    usetheseEarly = A.BlockPosition <= N & A.BlockPosition > 0;
    usetheseLate = A.BlockPosition >= -N & A.BlockPosition < 0;

    % Comparing rat and model behavior
    if strcmp(ratList{rr}, exampleRatMdl)
        [ratBlk, mdlBlk] =...
            compare_rat_and_model(BestFitLate.(ratList{rr}));
    end

    % Block transition dymnamics
    [ltom(rr,:), htom(rr,:)] =...
        block_dynamics_latency(A, twin, smoothfactor, true, [], [], true);

    % Previous reward regression for early and late trials
    [betasEarly(rr,:), tauEarly(rr), ~, ~,...
        ~, expParamsEarly(rr,:)] =...
        regress_latency_vs_rew(A, nback, false, true, true, 0,...
        usetheseEarly);
    [betasLate(rr,:), tauLate(rr), ~, ~,...
        ~, expParamsLate(rr,:)] =...
        regress_latency_vs_rew(A, nback, false, true, true, 0,...
        usetheseLate);

    % Simulated regressions
    [AMdlEarly, AMdlLate] = deal(BestFitEarly.(ratList{rr}).All.ratTrial);
    AMdlEarly.ITI = BestFitEarly.(ratList{rr}).All.LatMdl;
    AMdlLate.ITI = BestFitLate.(ratList{rr}).All.LatMdl;

    betasEarlySim(rr,:) =...
        regress_latency_vs_rew(AMdlEarly, nback, false, true, true, 0);
    betasLateSim(rr,:) =...
        regress_latency_vs_rew(AMdlLate, nback, false, true, true, 0);

    % Model comparison
    paramsEarly = BestFitEarly.final_params_mean(rr,:);
    paramsLate = BestFitLate.final_params_mean(rr,:);

    [~, nLL_early_early_trial] =...
        lognNLL_latency_secondhalf(paramsEarly, A.ITI, A,...
        'VanillaAlpha', [1.5 0 N 0 90]);
    nLLEarlyEarly(rr) = mean(nLL_early_early_trial, 'omitnan');

    [~, nLL_early_late_trial] =...
        lognNLL_latency_secondhalf(paramsEarly, A.ITI, A,...
        'VanillaAlpha', [1.5 -N 0 0 90]);
    nLLEarlyLate(rr) = mean(nLL_early_late_trial, 'omitnan');

    [~, nLL_late_early_trial] =...
        lognNLL_latency_secondhalf(paramsLate, A.ITI, A,...
        'VanillaAlpha', [1.5 0 N 0 90]);
    nLLLateEarly(rr) = mean(nLL_late_early_trial, 'omitnan');

    [~, nLL_late_late_trial] =...
        lognNLL_latency_secondhalf(paramsLate, A.ITI, A,...
        'VanillaAlpha', [1.5 -N 0 0 90]);
    nLLLateLate(rr) = mean(nLL_late_late_trial, 'omitnan');
end

% Pull recovered learning rates
alphaEarly = BestFitEarly.final_params_mean(:,1);
alphaLate = BestFitLate.final_params_mean(:,1);
Nrats = sqrt(length(ratList)); % For SEM

%% Model simulations
rng(1114) % Set random seed for simulation
alpha = 0.1; % Set learning rate

% Generate reward vector
R = [randsample(3:5, 15, true),...
    randsample(1:5, 40, true),...
    randsample(1:3, 15, true)];

% Preallocate value 
V = nan(size(R));
V(1) = 4.5;

% Loop over trials
for tt = 1:length(R)-1
    V(tt+1) = V(tt) + alpha*(R(tt) - V(tt)); % update value via delta-rule
end

V = (V-min(V))./(max(V) - min(V)); % min-max normalize for visualization

%% Save Data
save([savedir, 'Figure2_Data'], 'ratList', 'twin', 'nback',...
    'exampleRatMdl', 'ltom', 'htom', 'betasEarly', 'betasLate',...
    'betasEarlySim', 'betasLateSim',...
    'tauEarly', 'tauLate', 'nLLEarlyEarly', 'nLLEarlyLate',...
    'nLLLateEarly', 'nLLLateLate', 'expParamsEarly', 'expParamsLate',...
    'ratBlk', 'mdlBlk', 'alphaEarly', 'alphaLate', 'Nrats', 'V')

