function processMuscimolData(muscimolDataPath, codePath, savePath)
% Process muscimol data. Must run for plotFigure2_muscimol, plotFigure_S4  

%INPUTS: 
%   muscimolPath = local path to muscimol behavior data downloaded from Zenodo
%   codePath = local path to where code was saved
%   savePath = local path to where you want to save the outputs of this function

%% Set up paths
s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codePath, 'IgnoreCase', ispc);

if ~onPath % only add code dir to path if it already isn't
    addpath(genpath(codePath))
end

%% names of muscimol behavior data files
ratListC = {'ratTrial_S001Control'; 'ratTrial_S002Control'; ...
    'ratTrial_S003Control'; 'ratTrial_S009Control'; 'ratTrial_S014Control'; ...
    'ratTrial_S017Control'; 'ratTrial_S018Control'; 'ratTrial_S027Control'; ...
    'ratTrial_S028Control'};

ratListM = {'ratTrial_S001Muscimol'; 'ratTrial_S002Muscimol'; ...
    'ratTrial_S003Muscimol'; 'ratTrial_S009Muscimol'; ...
    'ratTrial_S014Muscimol'; 'ratTrial_S017Muscimol'; 'ratTrial_S018Muscimol'; ...
    'ratTrial_S027Muscimol'; 'ratTrial_S028Muscimol'};

%% run analysis functions
twin = 40; %window used for block transition dynamics, blocks are 40 successful trials
binSize = 5; %average block transition dynamics over bins of 5 trials
smoothfactor = 5;
nrats = length(ratListC);
nback = 10; %number of previous trials to regress against for regress_wt_vs_rew

for rr = 1:nrats
    
    % load control and muscimol data files for each rat
    S1 = load(strcat(muscimolDataPath, ratListC{rr}, '.mat'));
    S2 = load(strcat(muscimolDataPath, ratListM{rr}, '.mat'));

    %detrend wait times to account for any satiation effects
    S1.A = detrendwt_SS(S1.A);
    S2.A = detrendwt_SS(S2.A);
    
    %wait time curves
    [hiC(rr), loC(rr), mixC(rr)] = wtcurves_SS(S1.A);
    [hiM(rr), loM(rr), mixM(rr)] = wtcurves_SS(S2.A);
    
    %wait time dynamics
    [control.ltom(rr,:), control.htom(rr,:), control.mtol(rr,:), ...
        control.mtoh(rr,:), ~, ~, ~] =...
        block_dynamics_wt_binTrials(S1.A, twin, binSize, smoothfactor);

    [muscimol.ltom(rr,:), muscimol.htom(rr,:), muscimol.mtol(rr,:), ...
        muscimol.mtoh(rr,:), ~, ~, ~] =...
        block_dynamics_wt_binTrials(S2.A, twin, binSize, smoothfactor);
    
    %Regress wait times vs reward to get slopes and offsets
    control.slopes.mix(rr,:) = regressVolume(S1.A, []);
    muscimol.slopes.mix(rr,:) = regressVolume(S2.A, []);

    control.slopes.lo(rr,:) = regressVolume(S1.A, [], 'low');
    muscimol.slopes.lo(rr,:) = regressVolume(S2.A, [], 'low');

    control.slopes.hi(rr,:) = regressVolume(S1.A, [], 'high');
    muscimol.slopes.hi(rr,:) = regressVolume(S2.A, [], 'high');

    %regress wait time vs previous rewards (Mah 2023)
    [control.regress(rr,:), ~] = regress_wt_vs_rew(S1.A, nback);
    [muscimol.regress(rr,:), ~] = regress_wt_vs_rew(S2.A, nback);

    %ITI
    [hi_iti_C(rr), lo_iti_C(rr), mix_iti_C(rr), all_iti_C(rr)] = ...
        iticurves_SS(S1.A, 0); %post-unrewarded trials (Mah 2023)
    [hi_iti_M(rr), lo_iti_M(rr), mix_iti_M(rr), all_iti_M(rr)] = ...
        iticurves_SS(S2.A, 0);

    %split post low and post high mixed blocks into quartiles
    S1.A.trial_num = cell2mat(arrayfun(@(x) 1:S1.A.ntrials(x), ...
        1:length(S1.A.ntrials), 'uniformoutput', false))'; %add trial_num field to run quartileAnalysis
    S1.A.trainingstage = repmat(9, length(S1.A.block), 1);

    S2.A.trial_num = cell2mat(arrayfun(@(x) 1:S2.A.ntrials(x), ...
        1:length(S2.A.ntrials), 'uniformoutput', false))';
    S2.A.trainingstage = repmat(9, length(S2.A.block), 1);

    [~, ~, control.postLow_q1(rr,:), control.postHigh_q1(rr, :)] = quartileAnalysis_SS(S1.A);

    [~, ~, muscimol.postLow_q1(rr,:), muscimol.postHigh_q1(rr, :)] = quartileAnalysis_SS(S2.A);

end

%wait time curves
control.hi = cell2mat(arrayfun(@(x) hiC(x).wt, 1:nrats, 'uniformoutput', false)');
control.lo = cell2mat(arrayfun(@(x) loC(x).wt, 1:nrats, 'uniformoutput', false)');
control.mix = cell2mat(arrayfun(@(x) mixC(x).wt, 1:nrats, 'uniformoutput', false)'); 

muscimol.hi = cell2mat(arrayfun(@(x) hiM(x).wt, 1:nrats, 'uniformoutput', false)');
muscimol.lo = cell2mat(arrayfun(@(x) loM(x).wt, 1:nrats, 'uniformoutput', false)');
muscimol.mix =  cell2mat(arrayfun(@(x) mixM(x).wt, 1:nrats, 'uniformoutput', false)');

%iti curves
control.hi_iti.z = cell2mat(arrayfun(@(x) hi_iti_C(x).z(1), 1:nrats, ...
    'UniformOutput', false)');
control.lo_iti.z = cell2mat(arrayfun(@(x) lo_iti_C(x).z(1), 1:nrats, ...
    'UniformOutput', false)');
control.mix_iti.z = cell2mat(arrayfun(@(x) mix_iti_C(x).z(1), 1:nrats, ...
    'UniformOutput', false)');
control.all_iti.z = cell2mat(arrayfun(@(x) all_iti_C(x).z(1), 1:nrats, ...
    'UniformOutput', false)');

control.hi_iti.raw = cell2mat(arrayfun(@(x) hi_iti_C(x).raw(1), 1:nrats, ...
    'UniformOutput', false)');
control.lo_iti.raw = cell2mat(arrayfun(@(x) lo_iti_C(x).raw(1), 1:nrats, ...
    'UniformOutput', false)');
control.mix_iti.raw = cell2mat(arrayfun(@(x) mix_iti_C(x).raw(1), 1:nrats, ...
    'UniformOutput', false)');
control.all_iti.raw = cell2mat(arrayfun(@(x) all_iti_C(x).raw(1), 1:nrats, ...
    'UniformOutput', false)');

muscimol.hi_iti.z = cell2mat(arrayfun(@(x) hi_iti_M(x).z(1), 1:nrats, ...
    'UniformOutput', false)');
muscimol.lo_iti.z = cell2mat(arrayfun(@(x) lo_iti_M(x).z(1), 1:nrats, ...
    'UniformOutput', false)');
muscimol.mix_iti.z = cell2mat(arrayfun(@(x) mix_iti_M(x).z(1), 1:nrats, ...
    'UniformOutput', false)');
muscimol.all_iti.z = cell2mat(arrayfun(@(x) all_iti_M(x).z(1), 1:nrats, ...
    'UniformOutput', false)');

muscimol.hi_iti.raw = cell2mat(arrayfun(@(x) hi_iti_M(x).raw(1), 1:nrats, ...
    'UniformOutput', false)');
muscimol.lo_iti.raw = cell2mat(arrayfun(@(x) lo_iti_M(x).raw(1), 1:nrats, ...
    'UniformOutput', false)');
muscimol.mix_iti.raw = cell2mat(arrayfun(@(x) mix_iti_M(x).raw(1), 1:nrats, ...
    'UniformOutput', false)');
muscimol.all_iti.raw = cell2mat(arrayfun(@(x) all_iti_M(x).raw(1), 1:nrats, ...
    'UniformOutput', false)');

%% save
save([savePath 'control.mat'], 'control')
save([savePath 'muscimol.mat'], 'muscimol')

