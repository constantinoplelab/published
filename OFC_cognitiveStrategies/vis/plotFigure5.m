function plotFigure5(ephysPath, incongruentPath_expert, incongruentPath_naive)
%Plots data for figure 5. 
% YOU MUST RUN THE FUNCTION processIncongruentResponseData FOR BOTH EXPERT
% AND NAIVE DATA TO RUN THIS FUNCTION

%INPUTS: 
%   ephysPath = path to ephys data. Path should contain folders for both
%       expert and naive data
%   incongruentPath_expert = path to saved incongruent response file
%       experts
%   incongruentPath_naive = path to saved incongruent response file for
%       naive


%% general

expertPath = [ephysPath, filesep, 'Expert', filesep];
naivePath = [ephysPath, filesep, 'Naive', filesep];

timeBin = [-4 8]; % time window around each event saved in the SU structs
x = timeBin(1):0.05:timeBin(2);
wndw = [0 0.5]; %window for determining significance
w = arrayfun(@(y) find(x == wndw(y)), 1:2); %indices of window for computing significance

sem = @(x) std(x,'omitnan')./sqrt(sum(any(~isnan(x), 2)));

setDefaultFigProps

%% load data

tcaInfo_expert = load([expertPath, 'TCA', filesep, 'info_8comp.mat']);
tcaInfo_expert = tcaInfo_expert.info;

tcaInfo_naive = load([naivePath, 'TCA', filesep, 'info_8comp.mat']);
tcaInfo_naive = tcaInfo_naive.info;

clusts = unique(tcaInfo_expert.TCAclust);
nc = length(clusts);

tca_wndw = [-1 1];

%compute average responses aligned to task events to make heatmap
disp('This will take about 2 minutes to run')
avgResp_expert = averageResponse(expertPath, tcaInfo_expert, tca_wndw, 1);
avgResp_naive = averageResponse(naivePath, tcaInfo_naive, tca_wndw, 1);
ne = size(avgResp_expert, 2); % number of events for heatmap

%map cells to tca clusters
map_T = arrayfun(@(x) find(tcaInfo_expert.TCAclust == x), clusts, 'uniformoutput',...
        false); %find cells in each cluster, includes cluster "0" 
map_N = arrayfun(@(x) find(tcaInfo_naive.TCAclust == x), clusts, 'uniformoutput',...
        false); %find cells in each cluster, includes cluster "0" 
naiveToExpert = [1, 2, 6, 4, 3, 5, 7, 9, 8]; %qualitative mapping of naive to expert (+1 for indexing)


labels = {'Trial start' 'LED on' 'LED off' 'Reward' 'Opt out'};

%incongruent responses for experts
expert = load([incongruentPath_expert, 'incongruentResponses.mat']);

%incongruent responses for naive
naive = load([incongruentPath_naive, 'incongruentResponses.mat']);

%make combo table to pull cell responses by cluster
expertTable = join(expert.info, tcaInfo_expert);
naiveTable = join(naive.info, tcaInfo_naive);

%%
%expert figure
figure; hold on
tiledlayout(nc, ne+1, "TileSpacing", 'tight')

for c = 1:nc  
    for jj = 1:ne
        %plot cluster event-aligned heatmap
        byClust = avgResp_expert{jj}(map_T{c}, :);
        
        nexttile
        hold on
        imagesc(tca_wndw(1):0.05:tca_wndw(2),...
            linspace(1,length(byClust)),byClust,[-.5 1]);
        ylim([0 size(byClust, 1)]);
        yl = ylim;
        clim([-0.5 1.5]) %set same color axis for each plot

        if c == nc
            xlabel(labels{jj})
        else
            set(gca, 'xtick', [])
        end

        if jj == 1
            yticks([0 yl(2)])
        else
            set(gca, 'ytick', [])
        end
        set(gca,'TickDir','out');
                
    end
    %now plot cluster incongruent trial response (block-significant cells)
    clustPref = expert.incongPref{4}(expertTable.TCAclust == c-1, :);
    clustPref_Avg = mean(clustPref, 1, 'omitnan');
    clustPref_Sem = sem(clustPref);
    
    clustNonpref = expert.incongNonpref{4}(expertTable.TCAclust == c-1, :);
    clustNonpref_Avg = mean(clustNonpref, 1, 'omitnan');
    clustNonpref_Sem = sem(clustNonpref);
    
    %check significance
    nonVec = mean(clustNonpref(:, w(1):w(2)), 2, 'omitnan');
    nonVec = nonVec(~isnan(nonVec));

    prefVec = mean(clustPref(:, w(1):w(2)), 2, 'omitnan');
    prefVec = prefVec(~isnan(prefVec));

    [p_T, ~] = permute_test(nonVec, prefVec, 1000);

    ax = nexttile;
    shadedErrorBar(x, clustPref_Avg, clustPref_Sem, 'lineprops', ...
        {'color', [0.65 0.65 0.65], 'linewidth', 1}) 
    hold on
    shadedErrorBar(x, clustNonpref_Avg, clustNonpref_Sem, 'lineprops', ...
        {'color', 'k', 'linewidth', 1})
    xlim([-1 3])
    ylim([0 25])
    plot([0 0], [0 25], '--k', 'linewidth', 0.25)
    hold on
    set(gca, 'ytick', [])
    if c == nc
        set(gca, 'xtick', [-1, 0, 3])
        xlabel('Reward')
    else
        set(gca, 'xtick', [])
    end
    text(2, 23, string(p_T*9)) %bonferroni correct for 9 comparisons
    set(gca, 'TickDir', 'out'); box off;
    ax.YRuler.TickLabelGapOffset = 1;
    
end
clim([-0.5 1.5])
colorbar(ax, 'ticks', [-0.5, 0, 1.5])
sgtitle('Expert')
set(gcf,'renderer','painter')
set(gcf, 'units', 'centimeters', 'position', [10 10 10 12])

%%
%naive

figure; hold on
tiledlayout(nc, ne+1, "TileSpacing", 'tight')

for c = 1:nc  
    for jj = 1:ne
        %plot cluster event-aligned heatmap
        clust = naiveToExpert(c);
        byClust = avgResp_naive{jj}(map_N{clust}, :);
        
        nexttile
        hold on
        imagesc(tca_wndw(1):0.05:tca_wndw(2),...
            linspace(1,length(byClust)),byClust,[-.5 1]);
        ylim([0 size(byClust, 1)]);
        yl = ylim;
        clim([-0.5 1.5]) %set same color axis for each plot

        if c == nc
            xlabel(labels{jj})
        else
            set(gca, 'xtick', [])
        end

        if jj == 1
            yticks([0 yl(2)])
        else
            set(gca, 'ytick', [])
        end
        set(gca,'TickDir','out');
                
    end
    %now plot cluster incongruent trial response (block-significant cells)
    clustPref = naive.incongPref{4}(naiveTable.TCAclust == clust-1, :);
    clustPref_Avg = mean(clustPref, 1, 'omitnan');
    clustPref_Sem = sem(clustPref);
    
    clustNonpref = naive.incongNonpref{4}(naiveTable.TCAclust == clust-1, :);
    clustNonpref_Avg = mean(clustNonpref, 1, 'omitnan');
    clustNonpref_Sem = sem(clustNonpref);
    
    %check significance
    nonVec = mean(clustNonpref(:, w(1):w(2)), 2, 'omitnan');
    nonVec = nonVec(~isnan(nonVec));

    prefVec = mean(clustPref(:, w(1):w(2)), 2, 'omitnan');
    prefVec = prefVec(~isnan(prefVec));

    [p_N, ~] = permute_test(nonVec, prefVec, 1000);

    ax = nexttile;
    shadedErrorBar(x, clustPref_Avg, clustPref_Sem, 'lineprops', ...
        {'color', [0.65 0.65 0.65], 'linewidth', 1}) 
    hold on
    shadedErrorBar(x, clustNonpref_Avg, clustNonpref_Sem, 'lineprops', ...
        {'color', 'k', 'linewidth', 1})
    xlim([-1 3])
    ylim([0 25])
    plot([0 0], [0 25], '--k', 'linewidth', 0.25)
    hold on
    set(gca, 'ytick', [])
    if c == nc
        set(gca, 'xtick', [-1, 0, 3])
        xlabel('Reward')
    else
        set(gca, 'xtick', [])
    end
    text(2, 23, string(p_N*9)) %bonferroni correct for 9 comparisons
    set(gca, 'TickDir', 'out'); box off;
    ax.YRuler.TickLabelGapOffset = 1;
    
end
clim([-0.5 1.5])
colorbar(ax, 'ticks', [-0.5, 0, 1.5])
sgtitle('Naive')
set(gcf,'renderer','painter')
set(gcf, 'units', 'centimeters', 'position', [10 10 10 12])

