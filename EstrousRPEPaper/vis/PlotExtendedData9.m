function PlotExtendedData9(datadir, codedir)
%PlotExtendedData9 - Plots Extended Data 4. 
% INPUTS:
%   datadir - Local directory where ProcessData_ExtendedData9.mat was saved after running ProcessData_ExtendedData9
%   codedir - Local directory where code (e.g., published/EstrousRPE/) was saved

%% Set up paths
s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codedir, 'IgnoreCase', ispc);
if ~onPath % only add code dir to path if it isn't already
    addpath(genpath(codedir))
end

%% Load data
load([datadir, 'ProcessData_ExtendedData9'], ...
    'sortedPvsD_DOWN_NAccbackground', 'sortedEvsD_DOWN_NAccbackground', ...
    'sortedEvsD_UP_NAccbackground', 'dopamine_uptake_genes', ...
    'sortedPvsD', 'sortedDAgenes', 'sortedEvsD')

%% PLOT %%
figure;
set(gcf, 'Color', [1 1 1], 'Units', 'Inches', 'Position', [0, 0, 17, 9],...
    'PaperUnits', 'Inches', 'PaperSize', [9, 5], 'renderer','painters')

%--------------------------------------------------------------------------
%ED9a. Proestrus vs diestrus: downregulated genes compared to
% all identified proteins in NAcc as background, using STRING v12
%--------------------------------------------------------------------------
nexttile
barh(sortedPvsD_DOWN_NAccbackground.('strength')(end-4:end), 0.7, 'k'); hold on
yticks(1:5)
yticklabels(sortedPvsD_DOWN_NAccbackground.('termDescription')(end-4:end))
xlabel('Strength')
ylabel('Gene ontology')
ylim([0 6])
grid off; axis square; set(gca, 'TickDir', 'out'); box off
title('a Decreased in Proestrus vs. Diestrus')
subtitle('background is all identified proteins')

%--------------------------------------------------------------------------
%ED9b. Estrus vs. diestrus: downregulated proteins compared to
% all identified proteins in NAcc as background, using STRING v12
%--------------------------------------------------------------------------
nexttile
barh(sortedEvsD_DOWN_NAccbackground.('strength')(end-4:end), 0.7, 'k'); hold on
yticks(1:5)
yticklabels(sortedEvsD_DOWN_NAccbackground.('termDescription')(end-4:end))
xlabel('Strength')
ylabel('Gene ontology')
ylim([0 6])
grid off; axis square; set(gca, 'TickDir', 'out'); box off
title('b Decreased in Estrus vs. Diestrus')
subtitle('background is all identified proteins')

%--------------------------------------------------------------------------
%ED9c. Estrus vs diestrus: upregulated proteins compared to
% all identified proteins in NAcc as background, using STRING v12
%--------------------------------------------------------------------------
nexttile
barh(sortedEvsD_UP_NAccbackground.('strength')(end-4:end), 0.7, 'k'); hold on
yticks(1:5)
yticklabels(sortedEvsD_UP_NAccbackground.('termDescription')(end-4:end))
xlabel('Strength')
ylabel('Gene ontology')
ylim([0 6])
grid off; axis square; set(gca, 'TickDir', 'out'); box off
title('c Increased in Estrus vs. Diestrus')
subtitle('background is all identified proteins')

%--------------------------------------------------------------------------
%ED9d. Heatmap of dopamine-related terms, differential expression in
%our samples. Highlighting dopamine reuptake as term to change.
%--------------------------------------------------------------------------
nexttile
hAxes = gca;
imagesc(hAxes,sortedPvsD)
colormap(hAxes,flipud(pink))
colorbar
clim([-0.5 0.1])
gene_idx = find(ismember(sortedDAgenes, dopamine_uptake_genes));
for gn = 1:length(gene_idx)
    thisgene = gene_idx(gn);
    genename = sortedDAgenes{thisgene};
    text(1,thisgene,genename)
end
xticks([])
yticks([])
title('d Proestrus vs Diestrus')

nexttile
hAxes = gca;
imagesc(hAxes,sortedEvsD)
colormap(hAxes,flipud(pink))
colorbar
clim([-0.5 0.1])
gene_idx = find(ismember(sortedDAgenes, dopamine_uptake_genes));
for gn = 1:length(gene_idx)
    thisgene = gene_idx(gn);
    genename = sortedDAgenes{thisgene};
    text(1,thisgene,genename)
end
xticks([])
yticks([])
title('d Estrus vs Diestrus')

end



