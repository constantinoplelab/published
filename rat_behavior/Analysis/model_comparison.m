function [deltaCVnLL, deltaAIC, deltaBIC, pvals, ratList] =...
    model_comparison(mdl1_projectname, mdl1, mdl2_projectname, mdl2,...
    nll_arg, aic_arg, bic_arg, plotarg, varargin)
% Performs model comparison. All differenes are mdl2 - mdl1.
% Inputs:
%   mdl1_projectname - Project name for model 1
%   mdl1 - Model 1 type (e.g., Bayes)
%   mdl2_projectname - Project name for model 2
%   mdl2 - Model 2 type (e.g., Bayes)
%   nll_arg - Whether to calculate delta test nLL
%   aic_arg - Whether to calculate AIC
%   bic_arg - Whether to calculate BIC
%   plotarg - Whether to plot histograms
%   varargin - If you already have BestFits loaded, you can add them here
% Outputs:
%   deltaCVnLL - difference in cross-validated nLL
%   deltaAIC - difference in AIC = 2*n_params + 2*nLL
%   deltaBIC - deiffereince in BIC = n_params*log(n_trials) + 2*nLL
%   ratList - list of rats that you tested

if ispc
    mdl1_fitdir = ['Z:\BehavioralModel\' mdl1 'ModelFits\'...
        mdl1_projectname '\'];
    mdl2_fitdir = ['Z:\BehavioralModel\' mdl2 'ModelFits\'...
        mdl2_projectname '\'];
else
    mdl1_fitdir = ['/Volumes/server/BehavioralModel/' mdl1 'ModelFits/'...
        mdl1_projectname '/'];
    mdl2_fitdir = ['/Volumes/server/BehavioralModel/' mdl2 'ModelFits/'...
        mdl2_projectname '/'];
end

if isempty(varargin)
    A = load([mdl1_fitdir 'BestFit.mat']);
    BestFit_mdl1 = A.BestFit;
    
    B = load([mdl2_fitdir 'BestFit.mat']);
    BestFit_mdl2 = B.BestFit;
else
    BestFit_mdl1 = varargin{1};
    BestFit_mdl2 = varargin{2};
end

fieldlist1 = fields(BestFit_mdl1);
israt1 = structfun(@isstruct, BestFit_mdl1);
ratList1 = fieldlist1(israt1);

fieldlist2 = fields(BestFit_mdl2);
israt2 = structfun(@isstruct, BestFit_mdl2);
ratList2 = fieldlist2(israt2);

ratList = intersect(ratList1, ratList2);

pvals = nan(1, 3);

%% CrossVal nLL

if nll_arg
    TestNll_mdl1 = cellfun(@(x) BestFit_mdl1.(x).Test.nLL, ratList);
    TestNll_mdl2 = cellfun(@(x) BestFit_mdl2.(x).Test.nLL, ratList);
    
    deltaCVnLL = TestNll_mdl2 - TestNll_mdl1;
    p_nll = signrank(TestNll_mdl1, TestNll_mdl2);
else
    deltaCVnLL = nan;
    p_nll = nan;
end

pvals(1) = p_nll;

%% AIC
if aic_arg
    n_params_mdl1 = size(BestFit_mdl1.(ratList{1}).final_params, 2);
    n_params_mdl2 = size(BestFit_mdl2.(ratList{1}).final_params, 2);
    
    AIC_mdl1 = 2*n_params_mdl1 +...
        2*cellfun(@(x) BestFit_mdl1.(x).Test.nLL, ratList);
    AIC_mdl2 = 2*n_params_mdl2 +...
        2*cellfun(@(x) BestFit_mdl2.(x).Test.nLL, ratList);
    
    deltaAIC = AIC_mdl2 - AIC_mdl1;
    p_aic = signrank(AIC_mdl1, AIC_mdl2);
else
    deltaAIC = nan;
    p_aic = nan;
end

pvals(2) = p_aic;

%% BIC

if bic_arg
    n_params_mdl1 = size(BestFit_mdl1.(ratList{1}).final_params, 2);
    n_params_mdl2 = size(BestFit_mdl2.(ratList{1}).final_params, 2);
    
    ntrials = cellfun(@(r)...
        length(BestFit_mdl1.(r).Test.ratTrial.reward), ratList);
    
    BIC_mdl1 = log(ntrials)*n_params_mdl1 +...
        2*cellfun(@(x) BestFit_mdl1.(x).Test.nLL, ratList);
    BIC_mdl2 = log(ntrials)*n_params_mdl2 +...
        2*cellfun(@(x) BestFit_mdl2.(x).Test.nLL, ratList);
    
    deltaBIC = BIC_mdl2 - BIC_mdl1;
    p_bic = signrank(BIC_mdl1, BIC_mdl2);
else
    deltaBIC = nan;
    p_bic = nan;
end

pvals(3) = p_bic;

%%
if plotarg
    
    figure
    subplot(1, 3, 1)
    h1 = histogram(deltaCVnLL);
    if nll_arg
        xline(mean(deltaCVnLL, 'omitnan'), 'r--', linewidth=2, alpha=1)
    end
    
    xlim([-max(abs(h1.BinEdges)) max(abs(h1.BinEdges))])
    
    xlabel([mdl2 ' - ' mdl1], Interpreter='none')
    title({'Cross-Validated nLL', num2str(p_nll)})
    set(gca, Box='off')
    
    subplot(1, 3, 2)
    h2 = histogram(deltaAIC);
    xline(mean(deltaAIC, 'omitnan'), 'r--', linewidth=2, alpha=1)
    if aic_arg
        xlim([-max(abs(h2.BinEdges)) max(abs(h2.BinEdges))])
    end
    
    xlabel([mdl2 ' - ' mdl1], Interpreter='none')
    title({'AIC', num2str(p_aic)})
    set(gca, Box='off')
    
    subplot(1, 3, 3)
    h3 = histogram(deltaBIC);
    xline(mean(deltaBIC, 'omitnan'), 'r--', linewidth=2, alpha=1)
    if bic_arg
        xlim([-max(abs(h3.BinEdges)) max(abs(h3.BinEdges))])
    end
    
    xlabel([mdl2 ' - ' mdl1], Interpreter='none')
    title({'BIC', num2str(p_bic)})
    set(gca, Box='off')
    
    set(gcf, Position=[162 600 1519 345])
end
end