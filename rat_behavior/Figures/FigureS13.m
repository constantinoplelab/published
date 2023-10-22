function FigureS13(datadir, codedir)
%FigureS12 Males and females have comparable wait time ratios (top) and 
% trial initiation time ratios (bottom)
%   datadir = directory of dataset
%   codedir = director of code

%% Set up paths
s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codedir, 'IgnoreCase', ispc);

if ~onPath % only add code dir to path if it already isn't
    addpath(genpath(codedir))
end

%% Load Data

A = load([datadir filesep 'ModelFits' filesep...
    'RecoveringGenParams_Bayes' filesep...
    'RecoverGenParams.mat']);
genParams1 = A.gen_params;
recoveredParams1 = A.recovered_params;
paramNames1 = {'k_{mi}', 'k_{hi}', 'k_{lo}', 'D'};

B = load([datadir filesep 'ModelFits' filesep...
    'RecoveringGenParams_Bayes_BeliefState' filesep...
    'RecoverGenParams.mat']);

genParams2 = B.gen_params;
recoveredParams2 = B.recovered_params;
paramNames2 = {'k_{mi}', 'k_{hi}', 'k_{lo}', 'D'};

C = load([datadir filesep 'ModelFits' filesep...
    'RecoveringGenParams_BayesSubOpt_OnlyLambda_flat' filesep...
    'RecoverGenParams.mat']);

genParams3 = C.gen_params;
recoveredParams3 = C.recovered_params;
paramNames3 = {'k_{mi}', 'k_{hi}', 'k_{lo}', 'D', 'lambda'};

D = load([datadir filesep 'ModelFits' filesep...
    'RecoveringGenParams_AlphaPH' filesep...
    'RecoverGenParams.mat']);

genParams4 = D.gen_params;
recoveredParams4 = D.recovered_params;
paramNames4 = {'alpha', 'D'};

%% Plot

figure
for ii = 1:4
    subplot(4, 5, ii)
    plot(genParams1(:,ii), recoveredParams1(:,ii),...
        'ko', markerfacecolor='k', markersize=3)

    xlim([0 1])
    ylim([0 1])

    title(paramNames1{ii})
    axis square
    set(gca, box='off')
end

for ii = 1:4
    subplot(4, 5, ii+5)
    plot(genParams2(:,ii), recoveredParams2(:,ii),...
        'ko', markerfacecolor='k', markersize=3)

    xlim([0 1])
    ylim([0 1])

    title(paramNames2{ii})
    axis square
    set(gca, box='off')

end

for ii = 1:5
    subplot(4, 5, ii+10)
    plot(genParams3(:,ii), recoveredParams3(:,ii),...
        'ko', markerfacecolor='k', markersize=3)

    xlim([0 1])
    ylim([0 1])

    title(paramNames3{ii})
    axis square
    set(gca, box='off')

end

for ii = 1:2
    subplot(4, 5, ii+15)
    plot(genParams4(:,ii), recoveredParams4(:,ii),...
        'ko', markerfacecolor='k', markersize=4)

    xlim([0 1])
    ylim([0 1])

    title(paramNames4{ii})
    axis square
    set(gca, box='off')

end

set(gcf, Position=[894 76 787 869])
end