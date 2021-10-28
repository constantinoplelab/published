function [n_LR,n_winloss, n_prewinloss, n_safeRisky, p_LR,p_winloss, p_prewinloss, p_safeRisky ] = ofc_Ttests(A)
%ofc_Ttest: a wrapper for selectiveNeuronCalc to tests a neurons
%selectivity, for all neurons
%
%input: A. concatednated data, e.g. concatdata_ofc_start.mat
%
%
%ouput:
% number code for selectivity and p values from t tests
%   n_LR: +1 higher rate for left choice, -1 lower rate for left, 0
%   nonselective
%
%   n_winloss: +1 higher rate win. -1 lower rate in. 0 nonselective
%
%   n_prewinloss: +1 higher rate prev. win, -1 lower rate prev win, 0
%   nonselsectivte
%   
%   n_safeRisky: +1 higher rate safe choice, -1 lower rate safe choice, 0
%   nonselective
%
%   p_RL, etc. p values for each test

%Safe risky selective; Safe = 1, risky = -1. nonselective = 0
isRightSafe= @(A,j) A{j}.right_prob==1;
isLeftSafe= @(A,j) A{j}.left_prob==1;    
%bothsafe= @(A,j) isRightSafe(A,j) & isLeftSafe(A,j);
choseSafe = @(A,j) (isRightSafe(A,j) & A{j}.went_right==1) | ...
                (isLeftSafe(A,j) & A{j}.went_right==0);
mask1 = @(A,j) choseSafe(A,j) & A{j}.hits==1;
mask2 = @(A,j) ~choseSafe(A,j) & A{j}.hits==1;
[n_safeRisky,p_safeRisky] = selectiveNeuronCalc(A,mask1,mask2);

%left right selective. right = 1, left = -1, nonselective = 0
mask1 = @(A,j) A{j}.went_right == 0;
mask2 = @(A,j) A{j}.went_right == 1; 
[n_LR, p_LR] = selectiveNeuronCalc(A,mask1,mask2);

%check previous reward vs. no reward  selectivity. reward = 1, loss = -1,
%nonselective = 0
prev_hit_fcn = @(A,j)[nan; A{j}.hits(1:end-1)]; 
mask1 = @(A,j) prev_hit_fcn(A,j) ==1  & ~isnan(A{j}.hits);
mask2 = @(A,j) prev_hit_fcn(A,j) ==0  & ~isnan(A{j}.hits);
[n_prewinloss, p_prewinloss] = selectiveNeuronCalc(A,mask1,mask2);

%check current reward vs. no reward  selectivity. reward = 1, loss = -1,
%nonselective = 0
mask1 = @(A,j) A{j}.hits == 1 & ~isnan(A{j}.hits);
mask2 = @(A,j) A{j}.hits == 0 & ~isnan(A{j}.hits); 
[n_winloss, p_winloss] = selectiveNeuronCalc(A,mask1,mask2);



