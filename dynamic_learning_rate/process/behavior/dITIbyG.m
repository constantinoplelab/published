function [dITIbyG_small_mean, dITIbyG_large_mean,...
    dITIbyG_small_sem, dITIbyG_large_sem, p] = dITIbyG(ITI, RPE, A, israt)

% Calculate delta-belief G
[~, ~, ~, Belief] =...
    GenerateSynthData_Bayes(nan(1, 4), A,...
    'logn', false, nan);
belief = Belief(1, ~isnan(Belief(1,:)));
g = [nan, 1./(1 - abs(belief(2:end) - belief(1:end-1)))];
G = nan(size(Belief, 2), 1); G(~isnan(Belief(1,:))) = g;

dITI = [nan; diff(ITI)]; % Find trial by trial change in ITI

% Different parameters for rats vs. models
if israt
    postvios = [false; A.vios(1:end-1)];
    usethese = A.block==1 & postvios;
else
    usethese = A.block==1;
end

% Calculate bounds for small and large G
smallG = G < prctile(G, 50);
largeG = G >= prctile(G, 50);

% Find RPE bounds
negBnds = [prctile(RPE, 10), prctile(RPE, 20)];
negRPE = RPE > negBnds(1) & RPE < negBnds(2);

posBnds = [prctile(RPE, 80) prctile(RPE, 90)];
posRPE = RPE > posBnds(1) & RPE < posBnds(2);

% Average dITI by small G (Negative and Positive RPE)
dITIbyG_small_mean =...
    [mean(dITI(usethese & smallG & negRPE), 'omitnan'),...
    mean(dITI(usethese & smallG & posRPE), 'omitnan')];
dITIbyG_small_sem =...
    [sem(dITI(usethese & smallG & negRPE), 'omitnan'),...
    sem(dITI(usethese & smallG & posRPE), 'omitnan')];

% Average dITI by large G (Negative and Positive RPE)
dITIbyG_large_mean =...
    [mean(dITI(usethese & largeG & negRPE), 'omitnan'),...
    mean(dITI(usethese & largeG & posRPE), 'omitnan')];
dITIbyG_large_sem =...
    [sem(dITI(usethese & largeG & negRPE), 'omitnan'),...
    sem(dITI(usethese & largeG & posRPE), 'omitnan')];
end