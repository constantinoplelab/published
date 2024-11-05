function PlotFigure5(datadir, codedir)
%PlotFigure5 - Plots Figure 5. 
% INPUTS:
%   datadir - Local directory where ProcessData_Figure5.mat was saved after running ProcessData_Figure5
%   codedir - Local directory where code (e.g., published/EstrousRPE/) was saved

%% Set up paths
s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codedir, 'IgnoreCase', ispc);
if ~onPath % only add code dir to path if it already isn't
    addpath(genpath(codedir))
end

%% Load data
load([datadir, 'ProcessData_Figure5'], 'MS_pvd', 'log2FC_pvd', 'y_pvd',...
    'txt1_pvd', 'txt2_pvd', 'txt3_pvd', 'stagescompared_pvd',...
    'stageslfc_pvd', 'MS_evd', 'log2FC_evd', 'y_evd',...
    'txt1_evd', 'txt2_evd', 'txt3_evd', 'stagescompared_evd',...
    'stageslfc_evd', 'Stages_SERT', 'AreaByStage_SERT', 'pvals_SERT',...
    'Stages_DAT', 'AreaByStage_DAT', 'pvals_DAT', 'EMdata',...
    'FracSIGIC', 'FracSIGMem', 'groups', 'groupID');

%% PLOT %%
figure;
set(gcf, 'Color', [1 1 1], 'Units', 'Inches', 'Position', [0, 0, 17, 9],...
    'PaperUnits', 'Inches', 'PaperSize', [9, 5], 'renderer','painters')

%set generally used variables%% 
stagecolors = {'#EB8023', '#5E1AB7'}; %Estrus, Diestrus

%--------------------------------------------------------------------------
%5a. Volcano plots
%--------------------------------------------------------------------------
%proestrus/diestrus
nexttile
purplecolor = '#5E1AB7';
thresh = -log10(.05);
p = plot(log2FC_pvd(y_pvd<=thresh), y_pvd(y_pvd<=thresh), '.', 'Color', [.5 .5 .5]); hold on
p.Color(4) = 0.9;
plot(log2FC_pvd(y_pvd>thresh), y_pvd(y_pvd>thresh), '.k');
ix1 = find(cellfun(@(x)strcmp(x, txt1_pvd), MS_pvd.Gene));
ix2 = find(cellfun(@(x)strcmp(x, txt2_pvd), MS_pvd.Gene));
ix3 = find(cellfun(@(x)strcmp(x, txt3_pvd), MS_pvd.Gene));
txt1_pvd = 'DAT';
txt2_pvd = 'SERT';
plot(log2FC_pvd(ix1), y_pvd(ix1), '.', 'Color', purplecolor, 'MarkerSize', 15);
text(log2FC_pvd(ix1), y_pvd(ix1), txt1_pvd, 'Color', purplecolor);
plot(log2FC_pvd(ix2), y_pvd(ix2), '.', 'Color', purplecolor, 'MarkerSize', 15);
text(log2FC_pvd(ix2), y_pvd(ix2), txt2_pvd, 'Color', purplecolor);
plot(log2FC_pvd(ix3), y_pvd(ix3), '.', 'Color', purplecolor, 'MarkerSize', 15);
text(log2FC_pvd(ix3), y_pvd(ix3), txt3_pvd, 'Color', purplecolor);
title(['a ' stagescompared_pvd]);
xlabel(['Log-fold change [log2],',' ', stageslfc_pvd]);
ylabel('P-value [-log10]');
yticks(0:1:4)
xlim([-2.5 2.5]);
xticks(-2:1:2)
axis square; set(gca, 'TickDir', 'out'); box off

%estrus/diestrus
nexttile
thresh = -log10(.05);
p = plot(log2FC_evd(y_evd<=thresh), y_evd(y_evd<=thresh), '.', 'Color', [.5 .5 .5]); hold on
p.Color(4) = 0.9;
plot(log2FC_evd(y_evd>thresh), y_evd(y_evd>thresh), '.k');
ix1 = find(cellfun(@(x)strcmp(x, txt1_evd), MS_evd.Gene));
ix2 = find(cellfun(@(x)strcmp(x, txt2_evd), MS_evd.Gene));
ix3 = find(cellfun(@(x)strcmp(x, txt3_evd), MS_evd.Gene));
txt1_evd = 'SERT';
txt3_evd = 'DAT';
plot(log2FC_evd(ix1), y_evd(ix1), '.', 'Color', purplecolor, 'MarkerSize', 15);
text(log2FC_evd(ix1), y_evd(ix1), txt1_evd, 'Color', purplecolor);
plot(log2FC_evd(ix2), y_evd(ix2), '.', 'Color', purplecolor, 'MarkerSize', 15);
text(log2FC_evd(ix2), y_evd(ix2), txt2_evd, 'Color', purplecolor);
plot(log2FC_evd(ix3), y_evd(ix3), '.', 'Color', purplecolor, 'MarkerSize', 15);
text(log2FC_evd(ix3), y_evd(ix3), txt3_evd, 'Color', purplecolor);
title(['a ' stagescompared_evd]);
xlabel(['Log-fold change [log2],',' ', stageslfc_evd]);
ylabel('P-value [-log10]');
yticks(0:1:4)
xlim([-2.5 2.5]);
xticks(-2:1:2)
axis square; set(gca, 'TickDir', 'out'); box off

%--------------------------------------------------------------------------
%5c. SERT pixel area
%--------------------------------------------------------------------------
nexttile
for s=1:length(Stages_SERT)
    areadata = AreaByStage_SERT(:, s);
    avg = mean(areadata, 'omitnan');
    this_sem = std(areadata, 'omitnan')./sqrt(sum(~isnan(areadata)));
    errorbar(s, avg, this_sem, LineWidth=0.5,...
        color='k', capsize=10); hold on
    plot(s, avg, '.', MarkerSize=30,...
        color='k'); hold on
end
xticks([1 2 3])
xlim([0.5 3.5])
xticklabels(Stages_SERT)
ylabel('Area (pixels^2)')
yticks(0:0.4*10^4:1.8*10^4)
axis square; set(gca, 'TickDir', 'out'); box off
title(['c kruskal wallis p=' num2str(pvals_SERT.area.maineffect)])
subtitle(['posthoc PvsE p=' num2str(pvals_SERT.area.ProestrusEstrus),...
    ', PvsD p=' num2str(pvals_SERT.area.ProestrusDiestrus),...
    ', EvsD p=' num2str(pvals_SERT.area.EstrusDiestrus)])
effsize_pvd = effsize_ind(AreaByStage_SERT(:, 1), AreaByStage_SERT(:, 3));
effsize_evd = effsize_ind(AreaByStage_SERT(:, 2), AreaByStage_SERT(:, 3));
effsize_pve = effsize_ind(AreaByStage_SERT(:, 1), AreaByStage_SERT(:, 2));
disp(['SERT p vs d d=' num2str(effsize_pvd) ', e vs d d=' num2str(effsize_evd)...
    ', p vs e d=' num2str(effsize_pve)])

%--------------------------------------------------------------------------
%5d. DAT pixel area
%--------------------------------------------------------------------------
nexttile
for s=1:length(Stages_DAT)
    areadata = AreaByStage_DAT(:, s);
    avg = mean(areadata, 'omitnan');
    this_sem = std(areadata, 'omitnan')./sqrt(sum(~isnan(areadata)));
    errorbar(s, avg, this_sem, LineWidth=0.5,...
        color='k', capsize=10); hold on
    plot(s, avg, '.', MarkerSize=30,...
        color='k'); hold on
end
xticks([1 2])
xlim([0.5 2.5])
xticklabels(Stages_DAT)
yticks(0.8*10^4:0.4*10^4:2.2*10^4)
ylim([0.8*10^4 2.2*10^4])
ylabel('Area (pixels^2)')
axis square; set(gca, 'TickDir', 'out'); box off
title(['d rank sum p=' num2str(pvals_DAT.area.maineffect)])
effsize_pvd = effsize_ind(AreaByStage_DAT(:, 1), AreaByStage_DAT(:, 2));
disp(['DAT p vs d d=' num2str(effsize_pvd)])

%--------------------------------------------------------------------------
%5f. Membranous vs intracellular DAT particles identified with electron
% microscopy and their distance from the membrane (%membranous (includes
% ambiguous category because N-terminus of DAT is intracellular so could 
% still be bound to membrane if antibody appears near membrane)
%--------------------------------------------------------------------------
nexttile
FracSIGMem_bygroup = NaN(6, length(groups));
for s = 1:length(groups)
    groupidx = find(ismember(groupID, groups(s)));
    y = mean(FracSIGMem(groupidx, 1), 'omitnan')*100;
    this_sem = std(FracSIGMem(groupidx, 1), 'omitnan')./sqrt(size(FracSIGMem(groupidx, 1), 1))*100;
    ratsingroup = find(ismember(groupID, groups(s)));
    for rat = 1:length(ratsingroup)
        thisrat = ratsingroup(rat);
        plot(s, FracSIGMem(thisrat, 1)*100, '.', markersize=15, color=[0 0 0 0.3]); hold on
    end
    plot(s, y, '.k', markersize=30); hold on
    errorbar(s, y, this_sem, LineWidth=0.5, color='k', capsize=10); hold on
    ylabel('Percentage (%)')
    FracSIGMem_bygroup(:, s) = FracSIGMem(groupidx, 1);
end
pval = ranksum(FracSIGMem_bygroup(:, 1), FracSIGMem_bygroup(:, 2)); %only need to run this once because if not in one, then in the other
title(['f Membranous ranksum: p = ' num2str(pval)])
ylim([40 80])
yticks(40:10:80)
xticks(1:2)
xlim([0.5 2.5])
xticklabels({'Estrus', 'Diestrus'})
axis square; set(gca, 'TickDir', 'out'); box off

%intracellular
nexttile
for s = 1:length(groups)
    groupidx = find(ismember(groupID, groups(s)));
    y = mean(FracSIGIC(groupidx, 1), 'omitnan')*100;
    this_sem = std(FracSIGIC(groupidx, 1), 'omitnan')./sqrt(size(FracSIGIC(groupidx, 1), 1))*100;
    ratsingroup = find(ismember(groupID, groups(s)));
    for rat = 1:length(ratsingroup)
        thisrat = ratsingroup(rat);
        plot(s, FracSIGIC(thisrat, 1)*100, '.', markersize=15, color=[0 0 0 0.3]); hold on
    end
    plot(s, y, '.k', markersize=30); hold on
    errorbar(s, y, this_sem, LineWidth=0.5, color='k', capsize=10); hold on
    ylabel('Percentage (%)')
end
title('f Intracellular')
ylim([20 60])
yticks(20:10:60)
xticks(1:2)
xlim([0.5 2.5])
xticklabels({'Estrus', 'Diestrus'})
axis square; set(gca, 'TickDir', 'out'); box off

%CDF
%do not set ambiguous distnance to zero, only truly membranous SIGs were
%set to zero
nexttile
rawdata = EMdata.SIGMinimumDistancetoMembrane;
Estrus = rawdata(strcmp(EMdata.GroupID, 'estrus'));
Diestrus = rawdata(strcmp(EMdata.GroupID, 'diestrus'));
[~, kspval] = kstest2(Estrus, Diestrus);
h1 = cdfplot(Estrus); hold on
set(h1, 'Color', stagecolors{1});
set(h1, 'LineWidth', 0.75);
h2 = cdfplot(Diestrus); hold on
set(h2, 'Color', stagecolors{2});
set(h2, 'LineWidth', 0.75);
ylabel('Probability')
yticks(0:0.2:1)
title(['f Kolmogorov-Smirnov p = ' num2str(kspval)])
xlabel('SIG Distance from the Membrane')
axis square; set(gca, 'TickDir', 'out'); box off; grid off

%--------------------------------------------------------------------------
%5g-h. Simulation of DA fluorescence with reduced DAT
%--------------------------------------------------------------------------
simulate_DA_release
sgtitle('g-h')

end