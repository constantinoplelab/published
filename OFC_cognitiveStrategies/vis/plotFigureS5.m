function plotFigureS5(muscimolPath)

muscimolBehPath = [muscimolPath, filesep, 'BehaviorData', filesep];
muscimolPhysPath = [muscimolPath, filesep, 'Physiology', filesep];

ratListC = {'ratTrial_S001Control'; 'ratTrial_S002Control'; ...
    'ratTrial_S003Control'; 'ratTrial_S009Control'; 'ratTrial_S014Control'; ...
    'ratTrial_S017Control'; 'ratTrial_S018Control'; 'ratTrial_S027Control'; ...
    'ratTrial_S028Control'};

ratListM = {'ratTrial_S001Muscimol'; 'ratTrial_S002Muscimol'; ...
    'ratTrial_S003Muscimol'; 'ratTrial_S009Muscimol'; ...
    'ratTrial_S014Muscimol'; 'ratTrial_S017Muscimol'; 'ratTrial_S018Muscimol'; ...
    'ratTrial_S027Muscimol'; 'ratTrial_S028Muscimol'};

nrats = length(ratListC);

twin = 40;
binSize = 5;
smoothfactor = 5;

xvec = -twin:5:twin-5;

%general functions
sem = @(x) std(x,'omitnan')./sqrt(sum(any(~isnan(x), 2)));
setDefaultFigProps

%% Process muscimol rat data
[control, muscimol] = processMuscimolData(ratListC, ratListM, ...
    muscimolBehPath, twin, binSize, smoothfactor);

avgCm = mean(control.mix);
semCm = sem(control.mix);

avgEm = mean(muscimol.mix);
semEm = sem(muscimol.mix);

pCon = signrank(control.prevRew(:,1), control.prevRew(:,2)); 
pExp = signrank(muscimol.prevRew(:,1), muscimol.prevRew(:,2)); 

wt1 = control.hi(:,3)./control.lo(:,3);
wt2 = muscimol.hi(:,3)./muscimol.lo(:,3);

deltawt = wt2 - wt1;

bins = -0.5:0.025:0.5;

pslopes = signrank(control.slopes(:,2), muscimol.slopes(:,2));

%% process muscimol physiology data
out = muscimolPhysiology(muscimolPhysPath);

%% Regress WT vs reward

figure; 
tiledlayout(2, 4, 'TileSpacing', 'compact')

%example repressed cell during muscimol infusion
idx = find(out.xvec == -15):find(out.xvec == 30);

nexttile
plot(out.xvec(idx), out.hmat(idx), 'k', 'linewidth', 1.5)
hold on
plot([0 0], [0 10], '--k', 'linewidth', 1)
xlabel({'Time from infusion'; '(minutes)'})
ylabel('Firing rate (Hz)')
set(gca, 'TickDir', 'out'); box off;
ax1 = gca;
ax1.YRuler.TickLabelGapOffset = 1;
title('\rm Example cell')
axis square

%muscimol physiology summary
nexttile
plot(out.x, out.sigFit, 'color', [0.65 0.65 0.65], 'linewidth', 1.5)
hold on
scatter(out.bins, out.binnedratio, 100, '.k')
errorbar(out.bins, out.binnedratio, out.binnedEr, 'k', 'linestyle', 'none')
xl = get(gca, 'xlim');
plot([0 xl(2)], [1 1], '--k')
xlabel({'Distance from infusion'; 'site (mm)'})
ylabel('Normalized activity')
xlim([0 3.5])
ylim([0 2.5])
xticks(0:0.5:3.5)
text(0.2, 2.2, 'Sigmoid fit', 'color', [0.65 0.65 0.65], 'FontSize', 8)
set(gca, 'TickDir', 'out'); box off;
ax2 = gca;
ax2.YRuler.TickLabelGapOffset = 1;
axis square

% delta wt histogram

nexttile
h = histogram(deltawt, bins);
h.FaceColor = [0.65 0.65 0.65];
h.EdgeColor = [0.65 0.65 0.65];
xlim([-0.5 0.5])
ylim([0 3])
xline(0, 'k--', 'LineWidth', 1)
ylabel('N rats')    
xlabel({'\Delta wait time ratio'; '(muscimol - control)'})
axis square
set(gca, 'TickDir', 'out'); box off;
ax3 = gca;
ax3.YRuler.TickLabelGapOffset = 1;
box off   

%plot wait times in mixed blocks
x = 1:5;

nexttile
shadedErrorBar(x, avgCm, semCm, 'lineprops', {'color','k', ...
    'linewidth', 1.5})
shadedErrorBar(x, avgEm, semEm, 'lineprops', {'color',[0.65 0.65 0.65], ...
   'linewidth', 1.5})
set(gca, 'xTick', x);
set(gca, 'XTickLabels', {'5'; '10'; '20'; '40'; '80'});
set(gca, 'TickDir', 'out'); box off;
xlim([0 6]);
ylim([8 16])
xlabel('Reward offer')
ylabel('Mean wait time (s)')
text(0.2, 15.8, 'Muscimol', 'color', [0.65 0.65 0.65], 'FontSize', 8)
text(0.2, 14.5, 'Control', 'color', 'k', 'FontSize', 8)
text(3, 9, strcat('p =', num2str(pslopes)))
ax4 = gca;
ax4.YRuler.TickLabelGapOffset = 1;
axis square

nexttile
plot(0:6, control.regress(:, 1:7),...
    linewidth=0.5, color=[0.8 0.8 0.8])
shadedErrorBar(0:6, mean(control.regress(:,1:7)), sem(control.regress(:,1:7)),...
    lineprops={'color', 'k', 'linewidth', 2});
xticks(0:6)
xlim([-0.5, 6.5])
ylim([-0.3 0.6])
ax5 = gca;
ax5.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off;
axis square
xlabel('Trials back')
ylabel('Regression coefficient')
title('\rm Control')

nexttile
plot(0:6, muscimol.regress(:, 1:7),...
    linewidth=0.5, color=[0.8 0.8 0.8])
shadedErrorBar(0:6, mean(muscimol.regress(:,1:7)), sem(muscimol.regress(:,1:7)),...
    lineprops={'color', 'k', 'linewidth', 2});
xticks(0:6)
xlim([-0.5, 6.5])
ylim([-0.3 0.6])
ax6 = gca;
ax6.YRuler.TickLabelGapOffset = 1;
set(gca, 'TickDir', 'out'); box off;
axis square
xlabel('Trials back')
ylabel('Regression coefficient')
title('\rm Muscimol')

% conditional wait time for 20
x = [0.5 1];

nexttile
hold on
arrayfun(@(k) plot(x, control.prevRew(k,:) , 'color', [0 0 0 0.2]), 1:nrats)
plot(x, mean(control.prevRew, 'omitnan'), 'k', 'LineWidth', 2)
xlim([0.4 1.1])
ylim([-0.6 1])
xticks([0.5 1])
xticklabels({'<20' '>20'})
xlabel('Previous reward')
ylabel('Norm. wait time')
text(0.5, 0.75, strcat('p =', num2str(pCon)))
set(gca, 'TickDir', 'out'); box off;
ax7 = gca;
ax7.YRuler.TickLabelGapOffset = 1;
axis square
title('\rm Control')

nexttile
hold on
arrayfun(@(k) plot(x, muscimol.prevRew(k,:) , 'color', [0 0 0 0.2]), 1:nrats)
plot(x, mean(muscimol.prevRew, 'omitnan'), 'k', 'LineWidth', 2)
xlim([0.4 1.1])
ylim([-0.6 1])
xticks([0.5 1])
xticklabels({'<20' '>20'})
xlabel('Previous reward')
ylabel('Norm. wait time')
text(0.5, 0.75, strcat('p =', num2str(pExp)))
set(gca, 'TickDir', 'out'); box off;
ax8 = gca;
ax8.YRuler.TickLabelGapOffset = 1;
axis square
title('\rm Muscimol')
set(gcf,'renderer','painter')
set(gcf, 'units', 'centimeters', 'position', [10 10 20 10])
