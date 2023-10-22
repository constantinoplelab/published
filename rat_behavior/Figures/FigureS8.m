function FigureS8(datadir, codedir)
%FigureS8 -  Model comparison for wait times favors inferential over 
% retrospective model, but does not distinguish between inferential and 
% belief state models
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

%% Load data fits 
disp('loading data. may take a while.')
mdl1 = 'AlphaPH';
mdl1_projectname = 'FitAll_WT_AlphaPH_24July23_Final';

mdl2 = 'Bayes';
mdl2_projectname = 'FitAll_WT_Bayes_24July23_Final';

mdl4 = 'Bayes_BeliefState';
mdl4_projectname = 'FitAll_WT_BayesBeliefState_24July23_Final';

A = load([datadir 'ModelFits' filesep...
    mdl1_projectname filesep 'BestFit.mat']);
B = load([datadir 'ModelFits' filesep...
    mdl2_projectname filesep 'BestFit.mat']);
C = load([datadir 'ModelFits' filesep...
    mdl4_projectname filesep 'BestFit.mat']);

BestFit_mf = A.BestFit;
BestFit_mb = B.BestFit;
BestFit_mb_belief = C.BestFit;

%% Perform Model Comparison
[deltaNLL, deltaAIC, deltaBIC, pvals] =...
    model_comparison(mdl2_projectname, mdl2,...
    mdl1_projectname, mdl1,...
    true, true, true, false, BestFit_mb, BestFit_mf);

[deltaNLL2, deltaAIC2, deltaBIC2, pvals2] =...
    model_comparison(mdl2_projectname, mdl2,...
    mdl4_projectname, mdl4,...
    true, true, true, false, BestFit_mb, BestFit_mb_belief);

%% Plot histograms
figure
subplot(3, 2, 1)
h1 = histogram(deltaNLL);

xline(0, 'k--', linewidth=2, alpha=1)
xline(mean(deltaNLL, 'omitnan'), 'r-', linewidth=2, alpha=1)

set(h1, BinLimits = [-max(abs(h1.BinLimits)), max(abs(h1.BinLimits))])
set(gca, Box='off')

title({'Retro. - Infer.', ['nLL ' num2str(pvals(1))]})

subplot(3, 2, 3)
h2 = histogram(deltaAIC);

xline(0, 'k--', linewidth=2, alpha=1)
xline(mean(deltaAIC, 'omitnan'), 'r-', linewidth=2, alpha=1)

set(h2, BinLimits = [-max(abs(h2.BinLimits)), max(abs(h2.BinLimits))])
set(gca, Box='off')

title(['AIC ' num2str(pvals(2))])

subplot(3, 2, 5)
h2 = histogram(deltaBIC);

xline(0, 'k--', linewidth=2, alpha=1)
xline(mean(deltaBIC, 'omitnan'), 'r-', linewidth=2, alpha=1)

set(h2, BinLimits = [-max(abs(h2.BinLimits)), max(abs(h2.BinLimits))])
set(gca, Box='off')
title(['BIC ' num2str(pvals(3))])


subplot(3, 2, 2)
h1 = histogram(deltaNLL2);

xline(0, 'k--', linewidth=2, alpha=1)
xline(mean(deltaNLL2, 'omitnan'), 'r-', linewidth=2, alpha=1)

set(h1, BinLimits = [-max(abs(h1.BinLimits)), max(abs(h1.BinLimits))])
set(gca, Box='off')

title({'Belief. - Infer.', ['nLL ' num2str(pvals2(1))]})

subplot(3, 2, 4)
h2 = histogram(deltaAIC2);

xline(0, 'k--', linewidth=2, alpha=1)
xline(mean(deltaAIC2, 'omitnan'), 'r-', linewidth=2, alpha=1)

set(h2, BinLimits = [-max(abs(h2.BinLimits)), max(abs(h2.BinLimits))])
set(gca, Box='off')

title(['AIC ' num2str(pvals2(2))])

subplot(3, 2, 6)
h2 = histogram(deltaBIC2);

xline(0, 'k--', linewidth=2, alpha=1)
xline(mean(deltaBIC2, 'omitnan'), 'r-', linewidth=2, alpha=1)

set(h2, BinLimits = [-max(abs(h2.BinLimits)), max(abs(h2.BinLimits))])
set(gca, Box='off')
title(['BIC ' num2str(pvals2(3))])

set(gcf, position=[1000 140 681 805])
end