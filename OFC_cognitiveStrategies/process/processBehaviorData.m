function processBehaviorData(behaviorPath, codePath, savePath)
% Process behavior data. Must run for plotFigure1_expertBeh,
% plotFigure6_naiveBeh, plotFigureS1_rawWTexperts,
% plotFigureS9_naiveBehSupp

%INPUTS: 
%   behaviorPath = local path to behavior data downloaded from Zenodo
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

%% load data
load(strcat(behaviorPath, 'ratList.mat')); %list of rat data to use for behavior plots and model simulations
nrats = length(ratList);

behaviorPath_expert = behaviorPath;
behaviorPath_naive = [behaviorPath, filesep, 'All', filesep];


%% run analysis functions

disp('Running analysis functions')
disp('This will take ~5mins')

nsess = 15; %number of sessions considered "naive" 
twin = 40; %trial window for wait time dynamics plot
smoothfactor = 5; %window for causal smoothing
binSize = 1; % # of bins to average over
n = 0;

for rr = 1:nrats
    
    % load screened behavior data for expert rats
    E = load(strcat(behaviorPath_expert, 'ratTrial_', ratList{rr}, '.mat'));
    
    % load unscreened behavior data to pull out early training sessions
    N = load(strcat(behaviorPath_naive, 'ratTrialAll_', ratList{rr}, '.mat'));

    % use the first 15 sessions of data for early training/naive behavior
    N.A = splitData(N.A, nsess);

    % simulate divisive normalization agent
    divnorm_wt = divisiveNorm_fun_2(E.A, [50, 0.15, 60], 'log');
    div = E.A;
    div.wait_time = divnorm_wt;

    % simulate inference agent
    inf_wt = GenerateSynthData_Bayes_SS([0.25 0.30 0.2 .13], E.A, 'logn', 1, 8, 'log');
    inf = E.A;
    inf.wait_time = inf_wt;
    
    % detrend rat wait times to account for any satiation effects
    E.A = detrendwt_SS(E.A);
    N.A = detrendwt_SS(N.A);
    
    % wait time curves -- rat data
    [hiE(rr), loE(rr), mixE(rr)] = wtcurves_SS(E.A);
    [hiN(rr), loN(rr), mixN(rr)] = wtcurves_SS(N.A);
    
    % wait time dynamics -- rat data
    [expert.ltom(rr,:), expert.htom(rr,:), expert.mtol(rr,:), ...
        expert.mtoh(rr,:), expert.ltom_incong(rr,:), ...
        expert.htom_incong(rr,:), expert.dynamicsSems{rr}] =...
        block_dynamics_wt_binTrials(E.A, twin, binSize, smoothfactor);

    [naive.ltom(rr,:), naive.htom(rr,:), naive.mtol(rr,:), ...
        naive.mtoh(rr,:), naive.ltom_incong(rr,:), ...
        naive.htom_incong(rr,:), naive.dynamicsSems{rr}] =...
        block_dynamics_wt_binTrials(N.A, twin, binSize, smoothfactor);

    % wait time dynamics -- model data
    [infModel.ltom(rr,:), infModel.htom(rr,:), infModel.mtol(rr,:), ...
        infModel.mtoh(rr,:), ~, ~, ~] =...
        block_dynamics_wt_binTrials(inf, twin, binSize, smoothfactor);

    [divModel.ltom(rr,:), divModel.htom(rr,:), divModel.mtol(rr,:), ...
        divModel.mtoh(rr,:), ~, ~, ~] =...
        block_dynamics_wt_binTrials(div, twin, binSize, smoothfactor);
    
    % split post low and post high mixed blocks into quartiles
    [expert.postLow(rr,:), expert.postHigh(rr,:), expert.postLow_q1(rr,:),  ...
        expert.postHigh_q1(rr, :), expert.postLow_q1Con(rr,:), ...
        expert.postHigh_q1Con(rr,:)] = quartileAnalysis_SS(E.A);

    [naive.postLow(rr,:), naive.postHigh(rr,:), naive.postLow_q1(rr,:),  ...
        naive.postHigh_q1(rr, :), naive.postLow_q1Con(rr,:), ...
        naive.postHigh_q1Con(rr,:)] = quartileAnalysis_SS(N.A);

    [infModel.postLow(rr,:), infModel.postHigh(rr,:), infModel.postLow_q1(rr,:),  ...
        infModel.postHigh_q1(rr, :), infModel.postLow_q1Con(rr,:), ...
        infModel.postHigh_q1Con(rr,:)] = quartileAnalysis_SS(inf);

    [divModel.postLow(rr,:), divModel.postHigh(rr,:), divModel.postLow_q1(rr,:),  ...
        divModel.postHigh_q1(rr, :), divModel.postLow_q1Con(rr,:), ...
        divModel.postHigh_q1Con(rr,:)] = quartileAnalysis_SS(div);
    
    % Raw wait time transition dynamics
    % congruent volumes excluding 20
    [expert.ltomRaw(rr,:), expert.htomRaw(rr,:), expert.mtolRaw(rr,:), ...
        expert.mtohRaw(rr,:), ~, ~, ~] =...
        block_dynamics_wt_noNormalization(E.A, twin, binSize, 5, 1);

    [infModel.ltomRaw(rr,:), infModel.htomRaw(rr,:), infModel.mtolRaw(rr,:), ...
        infModel.mtohRaw(rr,:), ~, ~, ~] =...
        block_dynamics_wt_noNormalization(inf, twin, binSize, 5, 1);

    [naive.ltomRaw(rr,:), naive.htomRaw(rr,:), naive.mtolRaw(rr,:), ...
        naive.mtohRaw(rr,:), ~, ~, ~] =...
        block_dynamics_wt_noNormalization(N.A, twin, binSize, 5, 1);

    [divModel.ltomRaw(rr,:), divModel.htomRaw(rr,:), divModel.mtolRaw(rr,:), ...
        divModel.mtohRaw(rr,:), ~, ~, ~] =...
        block_dynamics_wt_noNormalization(div, twin, binSize, 5, 1);

    % 20ul only
    [expert.ltom20Raw(rr,:), expert.htom20Raw(rr,:), expert.mtol20Raw(rr,:), ...
        expert.mtoh20Raw(rr,:), ~, ~, ~] =...
        block_dynamics_wt_noNormalization(E.A, twin, binSize, 5, 0);

    [infModel.ltom20Raw(rr,:), infModel.htom20Raw(rr,:), infModel.mtol20Raw(rr,:), ...
        infModel.mtoh20Raw(rr,:), ~, ~, ~] =...
        block_dynamics_wt_noNormalization(inf, twin, binSize, 5, 0);

    [naive.ltom20Raw(rr,:), naive.htom20Raw(rr,:), naive.mtol20Raw(rr,:), ...
        naive.mtoh20Raw(rr,:), ~, ~, ~] =...
        block_dynamics_wt_noNormalization(N.A, twin, binSize, 2, 0);

    [divModel.ltom20Raw(rr,:), divModel.htom20Raw(rr,:), divModel.mtol20Raw(rr,:), ...
        divModel.mtoh20Raw(rr,:), ~, ~, ~] =...
        block_dynamics_wt_noNormalization(div, twin, binSize, 2, 0);

    % block duration examples
    examprats = [8 31];
    if rr == examprats(1) || rr == examprats(2)
        n = n+1;
        bdur = [];
        mintctr = 1;
        for k = 1:length(E.A.ntrials) %go through each session
            maxtctr = mintctr+E.A.ntrials(k)-1; %max trial index for the k session
            these = mintctr:maxtctr; %select trial indices for one session
            transitions = find(diff(E.A.block(these))~=0);
            bdur = [bdur; diff(transitions)];
            mintctr = mintctr+E.A.ntrials(k); %update counter for trial indices
        end
        expert.blockDuration{n} = bdur;
    end

end


%wait time curves
expert.hi = cell2mat(arrayfun(@(x) hiE(x).wt, 1:nrats, 'uniformoutput', false)');
expert.lo = cell2mat(arrayfun(@(x) loE(x).wt, 1:nrats, 'uniformoutput', false)');
expert.mix = cell2mat(arrayfun(@(x) mixE(x).wt, 1:nrats, 'uniformoutput', false)'); 

naive.hi = cell2mat(arrayfun(@(x) hiN(x).wt, 1:nrats, 'uniformoutput', false)');
naive.lo = cell2mat(arrayfun(@(x) loN(x).wt, 1:nrats, 'uniformoutput', false)');
naive.mix =  cell2mat(arrayfun(@(x) mixN(x).wt, 1:nrats, 'uniformoutput', false)');

save([savePath 'expert.mat'], 'expert')
save([savePath 'naive.mat'], 'naive')
save([savePath 'infModel.mat'], 'infModel')
save([savePath 'divModel.mat'], 'divModel')

