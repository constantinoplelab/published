function PlotFigure3(datadir, codedir)
%PlotFigure3 - Plots Figure 3. 
% INPUTS:
%   datadir - Local directory where ProcessData_Figure3.mat was saved after running ProcessData_Figure3
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
load([datadir, 'ProcessData_Figure3'],...
    'ratTrial', 'Kappa_trials', 'Lat_trials', 'betas_model',...
    'algn_vec', 'NAcc_ratlist', 'EncodingWindowAnalysis',...
    'pro_DA_binned', 'di_DA_binned', 'RPEbins', 'LatOpts', 'RPEs',...
    'RPEs_alphachange', 'delta', 'isblockchange', 'ltom', 'htom',...
    'mtol', 'mtoh', 'betas_model_3j', 'betas_RPE',...
    'betas_RPE_alphachange', 'ymean_RPEgain_binned', 'x_RPE_binned',...
    'BetasPro', 'BetasDi');

%% PLOT %%
figure;
set(gcf, 'Color', [1 1 1], 'Units', 'Inches', 'Position', [0, 0, 17, 9],...
    'PaperUnits', 'Inches', 'PaperSize', [9, 5], 'renderer','painters')   

%set generally used variables
cycle_colors = {'#E87003';'#7358C6'}; %orange (proestrus), purple (diestrus)
Stages = {'Proestrus', 'Diestrus'};
nback = 7;

%--------------------------------------------------------------------------
%3a. Reinforcement learning model that computes initiation times from kappa
%(expected value). Example predicted initiation time and kappa over blocks.
%--------------------------------------------------------------------------
subplot(6,3,1)
yl = [0 60];
mixedcolor=[0.4823529411764706 0.1568627450980392 0.48627450980392156];
fill([0 40 40 0], [yl(1) yl(1) yl(2) yl(2)],...
    mixedcolor, facealpha=0.15, edgecolor='none'); hold on
fill([40 80 80 40], [yl(1) yl(1) yl(2) yl(2)],...
    'r', facealpha=0.15, edgecolor='none'); hold on
fill([80 120 120 80], [yl(1) yl(1) yl(2) yl(2)],...
    mixedcolor, facealpha=0.15, edgecolor='none'); hold on
fill([120 160 160 120], [yl(1) yl(1) yl(2) yl(2)],...
    'b', facealpha=0.15, edgecolor='none'); hold on
fill([160 200 200 160], [yl(1) yl(1) yl(2) yl(2)],...
    mixedcolor, facealpha=0.15, edgecolor='none'); hold on
fill([200 240 240 200], [yl(1) yl(1) yl(2) yl(2)],...
    'r', facealpha=0.15, edgecolor='none'); hold on
fill([240 280 280 240], [yl(1) yl(1) yl(2) yl(2)],...
    mixedcolor, facealpha=0.15, edgecolor='none'); hold on
fill([280 320 320 280], [yl(1) yl(1) yl(2) yl(2)],...
    'b', facealpha=0.15, edgecolor='none'); hold on
fill([320 360 360 320], [yl(1) yl(1) yl(2) yl(2)],...
    mixedcolor, facealpha=0.15, edgecolor='none'); hold on
yyaxis left
plot(Kappa_trials(1:360), '-', color='k', linewidth=0.5); hold on
xlim([0 360])
ylim(yl)
yticks([])
xticks([])
ylabel('Reward rate')
yyaxis right
plot(Lat_trials(1:360), '-', color=[0 0 0 0.25], linewidth=0.5); hold on
yticks([])
xlabel('Trial')
ylabel('Modeled init. time')
set(gca, 'TickDir', 'out'); box off
title('a')

%--------------------------------------------------------------------------
%3b. Mean and SEM of predicted trial initiation times by block
%--------------------------------------------------------------------------
subplot(6,3,2)
iti_block_avg = NaN(1,2);
iti_block_er = NaN(1,2);
Block = ratTrial.block;
for bl = 2:3
    iti_block_avg(bl) = mean(LatOpt(Block==bl), 'omitnan');
    iti_block_er(bl) = std(LatOpt(Block==bl), 'omitnan')./...
        sqrt(sum(~isnan(LatOpt(Block==bl))));
end
y = [iti_block_avg(3) iti_block_avg(2)];
yer = [iti_block_er(3) iti_block_er(2)];
x = [1 2];
for bl = 1:2
    plot(x(bl), y(bl), '.k', 'MarkerSize', 30); hold on
    errorbar(x(bl), y(bl), yer(bl), color='k', LineWidth=0.5, capsize=0); hold on
end
xlim([.5 2.5]);
yticks([])
set(gca, 'TickDir', 'out'); box off
set(gca, 'xTick', x);
set(gca, 'xTickLabels', {'Low'; 'High'});
xlabel('Reward block');
ylabel('Modeled init. time');
set(gcf, 'Color', [1 1 1]);
set(gca, 'TickDir', 'out'); box off;
axis square
title('b')

%--------------------------------------------------------------------------
%3c. Regression of previous trial reward volumes on predicted initiation
%times
%--------------------------------------------------------------------------
subplot(6,3,3)
ratModel = ratTrial;
ratModel.ITI = LatOpt;
plot(1:nback, betas_model(1:nback), 'k', 'LineWidth', 0.5); hold on
xlim([0.5 6+0.5]) 
xticks(1:nback)
xticklabels(1:nback)
xlabel('Trials back')
ylabel('Regression coefficient')
ylim([-1.1 0.05])
yticks([])
yline(0, 'k--', 'linewidth', 0.5)
set(gca, 'TickDir', 'out'); box off;
axis square
title('c')

%--------------------------------------------------------------------------
%3d. Regression coefficient by time window between model-predicted RPE and 
% trial-by-trial DA z-score (last 20 trials) across time bins around every
% event. Use 100 ms bins. Pooling  across sessions for each rat, average 
% over rats. Includes 0.5 s baseline before CPIn.
%--------------------------------------------------------------------------
%Fit model to pooled sessions across rats and use that RPE
for aa  = 1:length(algn_vec)
    clear l
    % subplot(6,3,3+aa)
    subplot(1,length(algn_vec),aa)
    algn = algn_vec{aa};
    xtimes = EncodingWindowAnalysis.('G008').(algn).Times;
    ybetas = NaN(length(NAcc_ratlist), 150);
    for rat = 1:length(NAcc_ratlist)
        ybetas(rat, :) = EncodingWindowAnalysis.(NAcc_ratlist{rat}).(algn).Regression.RPE.betaMat(2,:); %first row is with the constant as the predictor
    end
    l(1) = shadedErrorBar(xtimes, mean(ybetas, 'omitnan'), sem(ybetas),...
        'lineProps', {'-', 'Color', 'k', 'LineWidth', 0.5}); %TRY WITH MEDIAN
    arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), l)
    xlabel('Time (s)')
    ylabel('Regression coefficient')
    axis square; grid off; set(gca, 'TickDir', 'out'); box off
    title(['d ' algn])
    ylim([-1 3.5])
    yticks(-1:1:3)
    xlim([-0.5 1.25])
    xline(0, '--k')
end

%--------------------------------------------------------------------------
%3e. DA AUC at offer cue as a function of model-predicted (using alpha fit to the 
% 10 trials) RPE
%--------------------------------------------------------------------------
subplot(6,3,10)
numbins = 7;
clear l
xvec = NaN(numbins, 1);
for b = 1:length(RPEbins)-1
    xvec(b) = (RPEbins(b) + RPEbins(b+1))./2;
end
l(1) = shadedErrorBar(xvec,...
    median(pro_DA_binned, 'omitnan'), sem(pro_DA_binned),...
    'lineProps', {'-', 'Color', cycle_colors{1}, 'LineWidth', 0.5}); hold on
l(1).mainLine.DisplayName = 'Proestrus';
l(2) = shadedErrorBar(xvec,...
    median(di_DA_binned, 'omitnan'), sem(di_DA_binned),...
    'lineProps', {'-', 'Color', cycle_colors{2}, 'LineWidth', 0.5}); hold on
l(2).mainLine.DisplayName = 'Diestrus';
legend('Location','best')
arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), l)
axis square; set(gca, 'TickDir', 'out'); box off
ylabel('DA AUC at Offer Cue')
yticks(-1:0.2:1)
xlabel('RPE (binned)')
xlim([xvec(1)-0.5 xvec(end)+0.5])
xticks(-3:1:3)
pvals = NaN(numbins-1, 1);
effsizes = NaN(numbins-1, 1);
for bin = 1:numbins
    pvals(bin, 1) = signrank(pro_DA_binned(:, bin), di_DA_binned(:, bin));
    if pvals(bin, 1) < 0.05
        text(xvec(bin), 0.65, '*')
    end
    effsizes(bin, 1) = effsize(pro_DA_binned(:, bin), di_DA_binned(:, bin));
end
title(['e ' num2str(numbins) ' bins, N=' num2str(length(NAcc_ratlist))])
disp(num2str(pvals'))
disp(num2str(effsizes'))

%--------------------------------------------------------------------------
%3f-j, l-m. Predicted block sensitivity and regression with added RPE to larger rewards
%--------------------------------------------------------------------------
num2plot = 200;
twin = 30;
xvec = (-twin:1:twin);
lowcolors = {'-b', '-k'};
highcolors = {'-r', '-k'};
for jj=1:2
    RPE = RPEs{jj};
    RPE_alphachange = RPEs_alphachange{jj};
    x = [1 2];

    %3g RPE over blocks of trials to show gain
    subplot(6,3,12)
    plot(1:num2plot, RPE(1:num2plot), color=cycle_colors{jj}); hold on
    ylabel('RPE'); xlabel('Trial')
    yticks([])
    set(gcf, 'Color', [1 1 1]);
    set(gca, 'TickDir', 'out'); box off;
    xline(isblockchange,...
        'k', 'LineWidth', 0.5); hold on
    xlim([0 num2plot])    
    set(gca, 'TickDir', 'out'); box off;    
    axis square
    title('g')

    %3h low - high initiation time
    subplot(6,3,13)
    plot(jj, delta{jj}, '.', color='k',...
        MarkerSize=30); hold on
    xlim([.5 2.5]);
    yticks([])
    ylabel('Low - high modeled init. time');
    set(gca, 'TickDir', 'out'); box off
    set(gca, 'xTick', x);
    set(gca, 'xTickLabels', {'Gain'; 'No gain'});
    set(gcf, 'Color', [1 1 1]);
    set(gca, 'TickDir', 'out'); box off;
    axis square
    title('h')

    % 3i latency dynamics
    subplot(6,3,15)
    % Plot as value transitions
    %group high to low value transitions
    hightolow = [htom{jj}; mtol{jj}]; %high to mix, mix to low
    hightolow_median = median(hightolow, 'omitnan');
    %group low to high value transitions
    %low to mix, mix to high
    lowtohigh = [ltom{jj}; mtoh{jj}];
    lowtohigh_median = median(lowtohigh, 'omitnan');
    shadedErrorBar(xvec, hightolow_median, ...
        sem(hightolow), 'lineprops', lowcolors{jj});
    shadedErrorBar(xvec, lowtohigh_median, ...
        sem(lowtohigh), 'lineprops', highcolors{jj});
    hold on
    ylims = ylim;
    line([0 0], [ylims(1) ylims(2)], 'Color', [0 0 0], 'LineStyle', '--');
    set(gca, 'TickDir', 'out'); box off
    xlabel('Trials from block switch');
    ylabel('Modeled init. time'); 
    yticks([])
    title('High to low value/Low to high value');
    xlim([-twin/3 twin])
    axis square
    set(gcf, 'Color', [1 1 1]);
    title('i')

    %3j initiation time regression coefficients
    subplot(6,3,14)
    ratModel = ratTrial;
    ratModel.ITI = LatOpts{jj};
    this_betas_model = betas_model_3j{jj};
    plot(1:nback, this_betas_model(1:nback), color=cycle_colors{jj}, LineWidth=0.5); hold on
    xlim([0.5 6.5]) 
    xticks(1:nback)
    xticklabels(1:nback)
    xlabel('Trials back')
    ylabel('Regression coefficient')
    yline(0, 'k--', 'linewidth', 0.5)
    set(gca, 'TickDir', 'out'); box off;
    axis square
    title('j Modeled init. time regression')

    %3l RPEs regression coefficients, RPE gain manipulation
    subplot(6,3,16)
    ratModelwRPE = ratTrial;
    ratModelwRPE.RPE = RPE;
    these_betas_RPE = betas_RPE{jj};
    plot(1:nback, these_betas_RPE(1:nback), color=cycle_colors{jj},...
        LineWidth=0.5); hold on
    xlim([0.5 nback + 0.5])
    xticks(1:nback)
    xticklabels(0:1:nback)
    xlabel('Trials back')
    ylabel('Regression coefficient')
    ylim([-17 23])
    yticks([])
    yline(0, 'k--', 'linewidth', 0.5)
    set(gca, 'TickDir', 'out'); box off;
    axis square
    title('l Modeled RPE regression, RPE gain manipulation')

    %3m RPEs regression coefficients, alpha manipulation
    subplot(6,3,17)
    ratModelwRPE = ratTrial;
    ratModelwRPE.RPE = RPE_alphachange;
    these_betas_RPE_alphachange = betas_RPE_alphachange{jj};
    plot(1:nback, these_betas_RPE_alphachange(1:nback),...
        color=cycle_colors{jj}, LineWidth=0.5); hold on
    xlim([0.5 nback + 0.5]) 
    xticks(1:nback)
    xticklabels(0:1:nback)
    xlabel('Trials back')
    ylabel('Regression coefficient')
    ylim([-14 19])
    yticks([])
    yline(0, 'k--', 'linewidth', 0.5)
    set(gca, 'TickDir', 'out'); box off;
    axis square
    title('m Modeled RPE regression, alpha manipulation')
    
end

%--------------------------------------------------------------------------
%3f Plot RPE as a function of gain applied to RPEs
%--------------------------------------------------------------------------
subplot(6,3,11)
plot(x_RPE_binned, ymean_RPEgain_binned,'.k',markersize=12); hold on
plot(min(RPEs{2}):1:max(RPEs{2}), min(RPEs{2}):1:max(RPEs{2}), '--k'); hold on
xlim([min(RPEs{2})+20 max(RPEs{2})+20])
ylim([min(RPEs{2})+20 max(RPEs{2})+20])
yticks([])
xticks([])
xlabel('RPE')
ylabel('RPE * \phi')
set(gca, 'TickDir', 'out'); box off;    
axis square   
title('f')

%--------------------------------------------------------------------------
%3k. Regression coefficients of DA as a function of current and previous 
% rewards in mixed blocks
%--------------------------------------------------------------------------
subplot(6,3,18)
clear l
l(1)=shadedErrorBar(1:nback, median(BetasPro(:, 2:nback+1), 'omitnan'),...
    sem(BetasPro(:, 2:nback+1)), 'lineProps', {'-', 'Color',...
    cycle_colors{1}, 'LineWidth', 1}); hold on
l(1).mainLine.DisplayName = Stages{1};
l(2)=shadedErrorBar(1:nback, median(BetasDi(:, 2:nback+1), 'omitnan'),...
    sem(BetasDi(:, 2:nback+1)), 'lineProps', {'-', 'Color',...
    cycle_colors{2}, 'LineWidth', 1}); hold on
l(2).mainLine.DisplayName = Stages{2};
legend('Location','best')
xlim([0.5 nback+0.5])
ylim([-0.07 0.255])
xticks(1:nback+1)
xticklabels(0:1:nback)
ylim([-0.06 0.29])
yline(0, '--', color = 'k', linewidth=0.75, HandleVisibility='off')
xlabel('Trials Back')
ylabel('Regression coefficient')
axis square; grid off; set(gca, 'TickDir', 'out'); box off
arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), l)
regresspvals = NaN(nback, 1);
effsizes = NaN(nback, 1);
for prevrew = 2:nback+1
    regresspvals(prevrew-1, 1) = signrank(BetasPro(:, prevrew),...
        BetasDi(:, prevrew));
    if regresspvals(prevrew-1, 1)<0.05
        text(prevrew-1, 0.25, '*', 'FontSize',14)
    end
    effsizes(prevrew-1, 1) = effsize(BetasPro(:, prevrew),...
        BetasDi(:, prevrew));
end
title(['k N=' num2str(length(NAcc_ratlist))])
disp(regresspvals)
disp(effsizes)

end