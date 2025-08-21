function PlotExtendedData5(datadir, codedir)
%PlotExtendedData5 - Plots Extended Data 5. 
% INPUTS:
%   datadir - Local directory where ProcessData_ExtendedData5.mat was saved after running ProcessData_ExtendedData5
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
load([datadir 'ProcessData_ExtendedData5'],...
    'NAcc_ratlist', 'Stages', 'AUC_byrat', 'pro_DA_binned',...
    'est_DA_binned', 'met_DA_binned', 'di_DA_binned',...
    'RewardDA_gfp', 'RewardDA_err_gfp', 'T', 'RewardDA_mCherry',...
    'RewardDA_err_mCherry', 'hi_gfp', 'lo_gfp', 'high_err_gfp',...
    'low_err_gfp','hi_mCherry', 'lo_mCherry', 'high_err_mCherry',...
    'low_err_mCherry', 'Alignments','rewards');

%% PLOT %%
figure;
set(gcf, 'Color', [1 1 1], 'Units', 'Inches', 'Position', [0, 0, 17, 9],...
    'PaperUnits', 'Inches', 'PaperSize', [9, 5], 'renderer','painters')

%--------------------------------------------------------------------------
%ED5a-b. Min-max normalized AUC for the largest reward (64 ul) during 
% mixed blocks over the cycle
%--------------------------------------------------------------------------
subplot(1, 3, 1)
AUCbystage = NaN(length(NAcc_ratlist), length(Stages));
for rat = 1:length(NAcc_ratlist)
    for s = 1:length(Stages)
        stagename = Stages{s};
        AUCbystage(rat, s) = AUC_byrat{rat}.(stagename).AUC_norm(1, end); hold on
    end
    plot(1:4, AUCbystage(rat, :), '-', color=[0 0 0 0.3])
end
xticks(1:1:4)
xlim([0.5 4.5])
xticklabels(Stages)
plot(1:4, mean(AUCbystage, 'omitnan'), '-k', linewidth=2)
thisplot = shadedErrorBar(1:4, mean(AUCbystage, 'omitnan'), sem(AUCbystage),...
    'lineprops', {'color', 'k', 'linewidth', 2});
set(thisplot.edge,'LineStyle','none')
arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), thisplot)
%stats
p_kw = kruskalwallis(AUCbystage,[],'off'); %comparison over the cycle
p_mvd = signrank(AUCbystage(:, 1), AUCbystage(:, 2));
p_mvp = signrank(AUCbystage(:, 1), AUCbystage(:, 3)); 
p_mve = signrank(AUCbystage(:, 1), AUCbystage(:, 4)); 
p_dvp = signrank(AUCbystage(:, 2), AUCbystage(:, 3));
p_dve = signrank(AUCbystage(:, 2), AUCbystage(:, 4));
p_pve = signrank(AUCbystage(:, 3), AUCbystage(:, 4));
subtitle({['kruskal wallis p=' num2str(p_kw)]...
    ['sign rank met vs di p=' num2str(p_mvd)]...
    ['sign rank met vs pro p=' num2str(p_mvp)]...
    ['sign rank met vs estrus p=' num2str(p_mve)]...
    ['sign rank di vs pro p=' num2str(p_dvp)]...
    ['sign rank di vs estrus p=' num2str(p_dve)]...
    ['sign rank pro vs estrus p=' num2str(p_pve)]...
    })
ylim([0.5 1.1])
yticks(0.4:0.2:1.2)
ylabel('Normalized DA AUC at offer cue')
set(gcf, 'Color', [1 1 1]);
set(gca, 'TickDir', 'out'); box off; axis square
title('a Mixed block 64 ul')
%--------------------------------------------------------------------------
%ED5a-b. Min-max normalized AUC for the largest reward (64 ul) during 
%high blocks over the cycle
%--------------------------------------------------------------------------
subplot(1, 3, 2)
AUCbystage = NaN(length(NAcc_ratlist), length(Stages));
for rat = 1:length(NAcc_ratlist)
    for s = 1:length(Stages)
        stagename = Stages{s};
        AUCbystage(rat, s) = AUC_byrat{rat}.(stagename).AUC_norm(2, end); hold on
    end
    plot(1:4, AUCbystage(rat, :), '-', color=[0 0 0 0.3])
end
xticks(1:1:4)
xlim([0.5 4.5])
xticklabels(Stages)
yticks(0.4:0.2:1.2)
plot(1:4, mean(AUCbystage, 'omitnan'), '-k', linewidth=2)
thisplot = shadedErrorBar(1:4, mean(AUCbystage, 'omitnan'), sem(AUCbystage),...
    'lineprops', {'color', 'k', 'linewidth', 2});
set(thisplot.edge,'LineStyle','none')
arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), thisplot)
%stats
p_kw = kruskalwallis(AUCbystage,[],'off'); %comparison over the cycle
p_mvd = signrank(AUCbystage(:, 1), AUCbystage(:, 2));
p_mvp = signrank(AUCbystage(:, 1), AUCbystage(:, 3)); 
p_mve = signrank(AUCbystage(:, 1), AUCbystage(:, 4)); 
p_dvp = signrank(AUCbystage(:, 2), AUCbystage(:, 3));
p_dve = signrank(AUCbystage(:, 2), AUCbystage(:, 4));
p_pve = signrank(AUCbystage(:, 3), AUCbystage(:, 4));
subtitle({['kruskal wallis p=' num2str(p_kw)]...
    ['sign rank met vs di p=' num2str(p_mvd)]...
    ['sign rank met vs pro p=' num2str(p_mvp)]...
    ['sign rank met vs estrus p=' num2str(p_mve)]...
    ['sign rank di vs pro p=' num2str(p_dvp)]...
    ['sign rank di vs estrus p=' num2str(p_dve)]...
    ['sign rank pro vs estrus p=' num2str(p_pve)]...
    })
set(gcf, 'Color', [1 1 1]);
set(gca, 'TickDir', 'out'); box off; axis square
title('b High block 64 ul')

%--------------------------------------------------------------------------
%ED5c. Min-max normalized AUC for the largest RPE over the cycle
%--------------------------------------------------------------------------
subplot(1, 3, 3)
met_DA_norm = met_DA_binned;
di_DA_norm = di_DA_binned;
pro_DA_norm = pro_DA_binned;
est_DA_norm = est_DA_binned;
for rat = 1:length(NAcc_ratlist)
   all_rat_means = [met_DA_binned(rat, :);...
   di_DA_binned(rat, :);...
   pro_DA_binned(rat, :);...
   est_DA_binned(rat, :)];
   max_for_rat = max(all_rat_means(:));
   min_for_rat = min(all_rat_means(:));

   normalized = (all_rat_means-min_for_rat)./(max_for_rat - min_for_rat);
   met_DA_norm(rat, :) = normalized(1, :);
   di_DA_norm(rat, :) = normalized(2, :);
   pro_DA_norm(rat, :) = normalized(3, :);
   est_DA_norm(rat, :) = normalized(4, :);
end
for rat = 1:length(NAcc_ratlist)
    plot(1:4, [met_DA_binned(rat, end) di_DA_norm(rat, end)...
        pro_DA_norm(rat, end) est_DA_norm(rat, end)],...
        '-', color=[0 0 0 0.3]); hold on
end
plot(1:4, [mean(met_DA_norm(:, end), 'omitnan')...
    mean(di_DA_norm(:, end), 'omitnan')...
    mean(pro_DA_norm(:, end), 'omitnan')...
    mean(est_DA_norm(:, end), 'omitnan')], '-k', linewidth=2)
thisplot = shadedErrorBar(1:4, [mean(met_DA_norm(:, end), 'omitnan')...
    mean(di_DA_norm(:, end), 'omitnan')...
    mean(pro_DA_norm(:, end), 'omitnan')...
    mean(est_DA_norm(:, end), 'omitnan')],...
    [sem(met_DA_norm(:, end))...
    sem(di_DA_norm(:, end))...
    sem(pro_DA_norm(:, end))...
    sem(est_DA_norm(:, end))],...
    'lineprops', {'color', 'k', 'linewidth', 2});
%stats
AUCbystage = [met_DA_norm(:, end), di_DA_norm(:, end),...
    pro_DA_norm(:, end), est_DA_norm(:, end)];
p_kw = kruskalwallis(AUCbystage,[],'off'); %comparison over the cycle
p_mvd = signrank(AUCbystage(:, 1), AUCbystage(:, 2));
p_mvp = signrank(AUCbystage(:, 1), AUCbystage(:, 3)); 
p_mve = signrank(AUCbystage(:, 1), AUCbystage(:, 4)); 
p_dvp = signrank(AUCbystage(:, 2), AUCbystage(:, 3));
p_dve = signrank(AUCbystage(:, 2), AUCbystage(:, 4));
p_pve = signrank(AUCbystage(:, 3), AUCbystage(:, 4));
subtitle({['kruskal wallis p=' num2str(p_kw)]...
    ['sign rank met vs di p=' num2str(p_mvd)]...
    ['sign rank met vs pro p=' num2str(p_mvp)]...
    ['sign rank met vs estrus p=' num2str(p_mve)]...
    ['sign rank di vs pro p=' num2str(p_dvp)]...
    ['sign rank di vs estrus p=' num2str(p_dve)]...
    ['sign rank pro vs estrus p=' num2str(p_pve)]...
    })
set(thisplot.edge,'LineStyle','none')
arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), thisplot)
xticks(1:1:4)
xlim([0.5 4.5])
xticklabels(Stages)
set(gcf, 'Color', [1 1 1]);
set(gca, 'TickDir', 'out'); box off; axis square
title('c Largest RPE bin')

%--------------------------------------------------------------------------
%ED5d. Motion-corrected GRABDA split by reward, aligned to each
% event
%--------------------------------------------------------------------------
figure;
set(gcf, 'Color', [1 1 1], 'Units', 'Inches', 'Position', [0, 0, 17, 9],...
    'PaperUnits', 'Inches', 'PaperSize', [9, 5], 'renderer','painters')
mycolors = {'#3852A3', '#4A3F98', '#7B287C', '#BD1F43', '#EC2024'};
rew_vols = {'4','8','16','32','64'};
for a = 1:length(Alignments)
    subplot(4, length(Alignments), a)
    hh = zeros(length(rewards), 1); %set up to name plots
    for rew = 1:length(rewards)
        da_mean = RewardDA_gfp{rew, a};
        standarderr = RewardDA_err_gfp{rew, a};
        err = [da_mean-standarderr...
            fliplr(da_mean+standarderr)];
        err(isnan(err)) = 0;
        hh(rew,1) = plot(T, da_mean, 'Color', mycolors{rew}, 'LineWidth', 0.5); hold on
        fl = fill([T fliplr(T)], err, 'k', 'FaceColor', mycolors{rew},...
            'LineStyle', 'none'); hold on
        set(fl, 'facealpha', 0.25);
        set(gca, 'TickDir', 'out'); box off
        xlabel(['Time from ', Alignments{a}]);
        xlim([-0.5 1.25])
        ylim([-0.5 1])
        axis square
        xline(0, '--k')
    end
    if a==1
        title('d Motion-corrected GRABDA')
    end
    if a == length(Alignments)
        legend(hh, rew_vols, location = 'best')
    end
end

%--------------------------------------------------------------------------
%ED5e. mCherry split by reward, aligned to each event
%--------------------------------------------------------------------------
for a = 1:length(Alignments)
    subplot(4, length(Alignments), length(Alignments)+a)
    for rew = 1:length(rewards)
        da_mean = RewardDA_mCherry{rew, a};
        standarderr = RewardDA_err_mCherry{rew, a};
        err = [da_mean-standarderr...
            fliplr(da_mean+standarderr)];
        err(isnan(err)) = 0;
        plot(T, da_mean, 'Color', mycolors{rew}, 'LineWidth', 0.5); hold on
        fl = fill([T fliplr(T)], err, 'k', 'FaceColor', mycolors{rew},...
            'LineStyle', 'none');
        set(fl, 'facealpha', 0.25);
        set(gca, 'TickDir', 'out'); box off
        xlabel(['Time from ', Alignments{a}]);
        xlim([-0.5 1.25])
        ylim([-0.5 1])
        axis square
        xline(0, '--k')
    end
    if a==1
        title('e mCherry')
    end
end

%--------------------------------------------------------------------------
%ED5f. Motion-corrected GRABDA split by block, aligned to each event
%--------------------------------------------------------------------------
cList = {'k'; 'r'; 'b'};
for a = 1:length(Alignments)
    subplot(4, length(Alignments), length(Alignments)*2+a)
    shadedErrorBar(T, hi_gfp{a}, high_err_gfp{a}, 'lineprops',...
        {'color', cList{2},'LineWidth', 0.5,'DisplayName', 'high'}); hold on
    shadedErrorBar(T, lo_gfp{a}, low_err_gfp{a}, 'lineprops',...
        {'color', cList{3},'LineWidth', 0.5,'DisplayName', 'low'}); hold on
    xlabel(['Time from ' Alignments{a} '(s)'])
    xline(0, '--k', 'HandleVisibility','off')
    set(gca,'TickDir','out'); box off
    ylim([-0.5 1])
    xlim([-0.5 1.25]);
    xticks(-1:0.5:1)
    axis square
    if a==1
        title('f Motion-corrected GRABDA')
    end
    if a == length(Alignments)
        legend('Location','best')
    end
end

%--------------------------------------------------------------------------
%ED5g. mCherry split by block, aligned to each event 
%--------------------------------------------------------------------------
for a = 1:length(Alignments)
    subplot(4, length(Alignments), length(Alignments)*3+a)
    shadedErrorBar(T, hi_mCherry{a}, high_err_mCherry{a}, 'lineprops',...
        {'color', cList{2},'LineWidth', 0.5,'DisplayName', 'high'}); hold on
    shadedErrorBar(T, lo_mCherry{a}, low_err_mCherry{a}, 'lineprops',...
        {'color', cList{3},'LineWidth', 0.5,'DisplayName', 'low'}); hold on
    xlabel(['Time from ' Alignments{a} '(s)'])
    xline(0, '--k')
    set(gca,'TickDir','out'); box off
    ylim([-0.5 1])
    xlim([-0.5 1.25]);
    xticks(-1:0.5:1)
    axis square
    if a==1
        title('g mCherry')
    end
end

end