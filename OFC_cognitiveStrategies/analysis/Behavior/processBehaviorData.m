function [expert, naive, infModel, divModel] = processBehaviorData(behaviorPath)

load(strcat(behaviorPath, 'ratList.mat')); %list of rat data to use for behavior plots and model simulations
nrats = length(ratList);

behaviorPath_expert = behaviorPath;
behaviorPath_naive = [behaviorPath, filesep, 'All', filesep];

nsess = 15; %number of sessions considered "naive" 
twin = 30; %trial window for wait time dynamics plot
smoothfactor = 5;
binSize = 1;

for rr = 1:nrats
    
    E = load(strcat(behaviorPath_expert, 'ratTrial_', ratList{rr}, '.mat'));
    N = load(strcat(behaviorPath_naive, 'ratTrialAll_', ratList{rr}, '.mat'));
    
    N.A = splitData(N.A, nsess); %naive data only includes the first nsess sessions

    % simulate divisive normalization agent
    divnorm_wt = divisiveNorm_fun(E.A, [5, 0.15, 60]);
    div = E.A;
    div.wait_time = divnorm_wt;

    % simulate inference agent
    inf_wt = GenerateSynthData_Bayes([0.23 0.3 0.2 .13 1], E.A, 'logn', 1, 8);
    inf = E.A;
    inf.wait_time = inf_wt;

    % wait time curves -- rat data
    [hiE(rr), loE(rr), mixE(rr)] = wtcurves(E.A);
    [hiN(rr), loN(rr), mixN(rr)] = wtcurves(N.A);
    
    % wait time dynamics -- rat data
    [expert.ltom(rr,:), expert.htom(rr,:), expert.mtol(rr,:), ...
        expert.mtoh(rr,:), ~, ~, ~] =...
        block_dynamics_wt_binTrials(E.A, twin, binSize, smoothfactor);

    [naive.ltom(rr,:), naive.htom(rr,:), naive.mtol(rr,:), ...
        naive.mtoh(rr,:), ~, ~, ~] =...
        block_dynamics_wt_binTrials(N.A, twin, binSize, smoothfactor);
    
    % wait time dynamics -- model data
    [infModel.ltom(rr,:), infModel.htom(rr,:), infModel.mtol(rr,:), ...
        infModel.mtoh(rr,:), ~, ~, ~] =...
        block_dynamics_wt_binTrials(inf, twin, binSize, smoothfactor);

    [divNorm.ltom(rr,:), divNorm.htom(rr,:), divNorm.mtol(rr,:), ...
        divNorm.mtoh(rr,:), ~, ~, ~] =...
        block_dynamics_wt_binTrials(div, twin, binSize, smoothfactor);

    % de-trend wait times
    E.A = detrendwt(E.A);
    N.A = detrendwt(N.A);
    
    %split post low and post high mixed blocks into quartiles
    [expert.postLow(rr,:), expert.postHigh(rr,:), expert.postLow_q1(rr,:),  ...
        expert.postHigh_q1(rr, :)] = quartileAnalysis(E.A);

    [naive.postLow(rr,:), naive.postHigh(rr,:), naive.postLow_q1(rr,:),  ...
        naive.postHigh_q1(rr, :)] = quartileAnalysis(N.A);

    [infModel.postLow(rr,:), infModel.postHigh(rr,:), infModel.postLow_q1(rr,:),  ...
        infModel.postHigh_q1(rr, :)] = quartileAnalysis(inf);

    [divModel.postLow(rr,:), divModel.postHigh(rr,:), divModel.postLow_q1(rr,:),  ...
        divModel.postHigh_q1(rr, :)] = quartileAnalysis(div);

end


%wait time curves
expert.hi = cell2mat(arrayfun(@(x) hiE(x).wt, 1:nrats, 'uniformoutput', false)');
expert.lo = cell2mat(arrayfun(@(x) loE(x).wt, 1:nrats, 'uniformoutput', false)');
expert.mix = cell2mat(arrayfun(@(x) mixE(x).wt, 1:nrats, 'uniformoutput', false)'); 

naive.hi = cell2mat(arrayfun(@(x) hiN(x).wt, 1:nrats, 'uniformoutput', false)');
naive.lo = cell2mat(arrayfun(@(x) loN(x).wt, 1:nrats, 'uniformoutput', false)');
naive.mix =  cell2mat(arrayfun(@(x) mixN(x).wt, 1:nrats, 'uniformoutput', false)');



