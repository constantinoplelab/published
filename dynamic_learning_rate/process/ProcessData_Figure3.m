function ProcessData_Figure3(BestFitEarly, codedir, savedir)
%ProcessData_Figure3 - function to process data to generate Figure 3. Saves
% data in savedir to be accessed by PlotFigure3
% INPUTS:
% BestFitEarly - BestFit.mat file found at: Mah_CellReports_BehavorialData/FitAll_ITI_VanillaAlpha_First10/BestFit.mat
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

fieldlist = fields(BestFitEarly);
ratList = fieldlist(structfun(@isstruct, BestFitEarly));

%% Parameters to save

iexample = 65; % example rat ID

%% Short example plot - for illustration
rng(22) % Set random seed

alpha0 = 0.3; % Learning rate
D = 20; % Scale

% Generate rewards
r = [randsample(3:5, 40, true),...
    randsample(1:5, 40, true),...
    randsample(1:3, 40, true)]';
r(41:42) = randsample(3:5, 2, true);
r(77:80) = randsample(3:5, 4, true);

% Pull model predicted gains
[~, ~, ~, ~, ~, G_DB1, BeliefExample] =...
    generate_deltabelief_mdl(r, alpha0, D);
[~, ~, ~, ~, ~, G_PH1] =...
    generate_ph_mdl(r, alpha0, D);
[~, ~, ~, ~, ~, G_MS1] =...
    generate_ms_mdl(r, alpha0, D);

%% Variance by block position

D = 7; % Scale parameter

% Preallocate data structures
[varByPosRat, varByPosPH, varByPosDB, varByPosMS] =...
    deal(nan(length(ratList), 2));

[itiRat, itiPH, itiDB, itiMS] = deal(cell(size(ratList)));

rng(724) % Set random seed for consistency
for rr = 1:length(ratList)
    fprintf('%d out of %d\n', rr, length(ratList)) % display progress

    % Pull rat data
    A = BestFitEarly.(ratList{rr}).All.ratTrial; 

    % threshold rat ITI
    itiRat{rr} = A.ITI;
    itiRat{rr}(A.trial_num == 1) = nan;
    itiRat{rr}(itiRat{rr} > prctile(itiRat{rr}, 90)) = nan;

    % Generate model ITIs
    alpha0 = rand*0.3; % Randomly generate learning rate

    R = convertreward(A.reward);
    itiPH{rr} =...
        generate_ph_mdl(R, alpha0, D);
    itiDB{rr} =...
        generate_deltabelief_mdl(R, alpha0, D);
    itiMS{rr} =...
        generate_ms_mdl(R, alpha0, D);

    % Logical vectors to pull early or late trials
    postvios = [false; A.vios(1:end-1)];
    usetheseEarly = A.block==1 & postvios &...
        A.BlockPosition > 0 & A.BlockPosition <= 10;
    usetheseLate = A.block==1 & postvios&...
        A.BlockPosition >= -10 & A.BlockPosition < 0;

    % Variances by early and late
    varByPosRat(rr,:) =...
        [var(itiRat{rr}(usetheseEarly), 'omitnan'),...
        var(itiRat{rr}(usetheseLate), 'omitnan')];
    varByPosPH(rr,:) =...
        [var(itiPH{rr}(usetheseEarly)), var(itiPH{rr}(usetheseLate))];
    varByPosDB(rr,:) =...
        [var(itiDB{rr}(usetheseEarly)), var(itiDB{rr}(usetheseLate))];
    varByPosMS(rr,:) =...
        [var(itiMS{rr}(usetheseEarly)), var(itiMS{rr}(usetheseLate))];
end

% Log-ratio for variance
varRatioRat = log(varByPosRat(:,2)./varByPosRat(:,1));
varRatioPH = log(varByPosPH(:,2)./varByPosPH(:,1));
varRatioDB = log(varByPosDB(:,2)./varByPosDB(:,1));
varRatioMS = log(varByPosMS(:,2)./varByPosMS(:,1));

%% Delta-ITI by delta-Belief

% preallocate
[dITIbyG_small, dITIbyG_large] = deal(nan(length(ratList), 2));
rpes = cell(size(ratList));

[negBnds, posBnds] = deal(nan(length(ratList), 2));

% Loop over rats
for rr = 1:length(ratList)
    fprintf('%d out of %d\n', rr, length(ratList)) % check progress
    A = BestFitEarly.(ratList{rr}).All.ratTrial; % pull data

    % Pull and filter data
    ITI = A.ITI;
    ITI(A.trial_num == 1) = nan; ITI(ITI > prctile(ITI, 99)) = nan;
    RPE = BestFitEarly.(ratList{rr}).All.RPE;
    RPE(abs(RPE) > 5) = nan;
    rpes{rr} = RPE;

    % Condition on delta-belief
    [dITIbyG_small(rr,:), dITIbyG_large(rr,:)] =...
        dITIbyG(ITI, rpes{rr}, A, true);
end

%% Model examples
% Rat data for model
A = BestFitEarly.(ratList{1}).All.ratTrial;
R = convertreward(A.reward);

% learning rate
alpha0 = 0.15;

% Generate data
[~, ITI_db, ~, RPE_db] = generate_deltabelief_mdl(R, alpha0, D);
[~, ITI_ph, ~, RPE_ph] = generate_ph_mdl(R, alpha0, D);
[~, ITI_ms, ~, RPE_ms] = generate_ms_mdl(R, alpha0, D);

% Condition on delta-belief
[dITIbyG_small_mean_DB, dITIbyG_large_mean_DB,...
    dITIbyG_small_sem_DB, dITIbyG_large_sem_DB] =...
    dITIbyG(ITI_db, RPE_db, A, false);

[dITIbyG_small_mean_PH, dITIbyG_large_mean_PH,...
    dITIbyG_small_sem_PH, dITIbyG_large_sem_PH] =...
    dITIbyG(ITI_ph, RPE_ph, A, false);

[dITIbyG_small_mean_MS, dITIbyG_large_mean_MS,...
    dITIbyG_small_sem_MS, dITIbyG_large_sem_MS] =...
    dITIbyG(ITI_ms, RPE_ms, A, false);

%% Changepoint model

rng(724) % Set random seed
R = generate_reward(50); % Generate rewards

% Calculate run-length posterior
RunLengthPosterior = calculate_changepoint_model(R, 75);

% Look at Bayes model performance on the same data
[~, ~, ~, ~, ~, ~, belief2compare] =...
    generate_deltabelief_mdl(...
    R, nan, nan);

% Probablility that a changepoint occured
change_point = RunLengthPosterior(2,:)>0.2;
change_point(2) = false;

% Get model gains
[~, ~, ~, ~, ~, G_db] = generate_deltabelief_mdl(R, alpha0, D);
[~, ~, ~, ~, ~, G_ph] = generate_ph_mdl(R, alpha0, D);
[~, ~, ~, ~, ~, G_ms] = generate_ms_mdl(R, alpha0, D);

% Look at gain conditioned by Changepoint
G_dbByChangePoint =...
    [mean(G_db(~change_point), 'omitnan'),...
    mean(G_db(change_point), 'omitnan') ;...
    sem(G_db(~change_point), 'omitnan'),...
    sem(G_db(change_point), 'omitnan')];

G_phByChangePoint =...
    [mean(G_ph(~change_point), 'omitnan'),...
    mean(G_ph(change_point), 'omitnan') ;...
    sem(G_ph(~change_point), 'omitnan'),...
    sem(G_ph(change_point), 'omitnan')];

G_msByChangePoint =...
    [mean(G_ms(~change_point), 'omitnan'),...
    mean(G_ms(change_point), 'omitnan') ;...
    sem(G_ms(~change_point), 'omitnan'),...
    sem(G_ms(change_point), 'omitnan')];

% Find BOCD block transitions (when max run-length decreases)
[~, i] = max(RunLengthPosterior);

BOCD_blockTransitions = find(diff(i) < 0);

% Find bayes model block transitions (When inferred block changes)
[~, blkInf] = max(belief2compare);
Bayes_blockTransitions = find(diff(blkInf) ~= 0);

% Find number of trials since last real block transitions
BOCD_nTrials =....
    arrayfun(@(n) min(abs(BOCD_blockTransitions - 40*n)), 1:199)';
Bayes_nTrials =....
    arrayfun(@(n) min(abs(Bayes_blockTransitions - 40*n)), 1:199)';

Bayes_nTrials_binned =...
    [arrayfun(@(n) mean(Bayes_nTrials(BOCD_nTrials==n)), 0:10);...
    arrayfun(@(n) std(Bayes_nTrials(BOCD_nTrials==n)), 0:10)];

% Correlate number of trials to infer block
[r, p] = corr(BOCD_nTrials, Bayes_nTrials);

%% Save Data

save([savedir, 'Figure3_Data'], 'BeliefExample',...
    'G_DB1', 'G_PH1', 'G_MS1', 'G_db', 'G_ph', 'G_ms',...
    'varRatioRat', 'varRatioPH', 'varRatioMS', 'varRatioDB',...
    'dITIbyG_small', 'dITIbyG_large', 'dITIbyG_small_mean_DB',...
    'dITIbyG_small_sem_DB', 'dITIbyG_large_mean_DB',...
    'dITIbyG_large_sem_DB', 'dITIbyG_small_mean_PH',...
    'dITIbyG_small_sem_PH', 'dITIbyG_large_mean_PH', ...
    'dITIbyG_large_sem_PH', 'dITIbyG_small_mean_MS',...
    'dITIbyG_large_mean_MS', 'dITIbyG_small_sem_MS',...
    'dITIbyG_large_sem_MS', 'rpes', 'iexample', 'negBnds', 'posBnds',...
    'G_msByChangePoint', 'BOCD_nTrials',...
    'Bayes_nTrials_binned', 'RunLengthPosterior', 'G_dbByChangePoint',...
    'G_phByChangePoint', 'r', 'p');

end
