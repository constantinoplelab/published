function plotFigure4_TCA(ephysPath, codePath, processedEphysPath)
% Plot figure 4. Process and plot heatmat of neural responses based on TCA model
% Process and plot responses to 20ul trials pre- and post-incongruent trial
% by TCA cluster

s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codePath, 'IgnoreCase', ispc);

if ~onPath % only add code dir to path if it already isn't
    addpath(genpath(codePath))
end

%% general

timeBin = [-4 8]; % time window around each event saved in the SU structs
x = timeBin(1):0.05:timeBin(2);
wndw = [-0.5 1];
w = arrayfun(@(y) find(x == wndw(y)), 1:2); %indices of window for computing significance
xvec = w(1):w(2);

sem = @(x) std(x,'omitnan')./sqrt(sum(any(~isnan(x), 2)));

setDefaultFigProps

%% OFC path
OFCpath = [ephysPath, filesep, 'Expert', filesep, 'OFC', filesep];

%% load incongruent responses
incongruentResponses = load([processedEphysPath, 'incongruentResponses.mat']);

%% load data

tcaInfo_expert = load([OFCpath, 'TCA', filesep, 'info_8comp.mat']);
tcaInfo_expert = tcaInfo_expert.info;

clusts = unique(tcaInfo_expert.TCAclust);
nc = length(clusts);

tca_wndw = [-1 1];

%compute average responses aligned to task events to make heatmap
disp('This will take about 2 minutes to run')
avgResp_expert = averageResponse(OFCpath, tcaInfo_expert, tca_wndw, 1);
ne = size(avgResp_expert, 2); % number of events for heatmap

%map cells to tca clusters
map_T = arrayfun(@(x) find(tcaInfo_expert.TCAclust == x), clusts, 'uniformoutput',...
        false); %find cells in each cluster, includes cluster "0" 

%make combo table to pull cell responses by cluster
expertTable = join(incongruentResponses.info, tcaInfo_expert);


%% plot

figure; hold on
tiledlayout(nc, ne+4, "TileSpacing", 'tight')

labels = {'Trial start' 'LED on' 'LED off' 'Reward', 'Opt out'};
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

    % now plot cluster incongruent responses for trial start and reward
    for  a = [1 4]

        clustPre_pref = incongruentResponses.pre_pref{a}(expertTable.TCAclust == c-1, :);
        clustPre_pref_Avg = mean(clustPre_pref, 1, 'omitnan');
        clustPre_pref_Sem = sem(clustPre_pref);

        clustPre_nonpref = incongruentResponses.pre_nonpref{a}(expertTable.TCAclust == c-1, :);
        clustPre_nonpref_Avg = mean(clustPre_nonpref, 1, 'omitnan');
        clustPre_nonpref_Sem = sem(clustPre_nonpref);

        clustPost_pref = incongruentResponses.post_pref{a}(expertTable.TCAclust == c-1, :);
        clustPost_pref_Avg = mean(clustPost_pref, 1, 'omitnan');
        clustPost_pref_Sem = sem(clustPost_pref);

        clustPost_nonpref = incongruentResponses.post_nonpref{a}(expertTable.TCAclust == c-1, :);
        clustPost_nonpref_Avg = mean(clustPost_nonpref, 1, 'omitnan');
        clustPost_nonpref_Sem = sem(clustPost_nonpref);

        % run permutation test on each time bin
        pExp_post = arrayfun(@(y) permute_test(clustPost_nonpref(:,y), ...
            clustPost_pref(:,y), 1000, 1), xvec);
        pExp_pre = arrayfun(@(y) permute_test(clustPre_nonpref(:,y), ...
            clustPre_pref(:,y), 1000, 1), xvec);

        % find significant bins
        pExp_pre_sig = pExp_pre < 0.05;
        pExp_post_sig = pExp_post < 0.05;


        %find 3 consecutive significant bins
        sig_pre = false(1, length(xvec));
        a0 = [0 pExp_pre_sig 0]; % add 0 at beginning and end to find pattern. also shifts index so start is the first true value
        startTrue = strfind(a0,[0 1]);
        endTrue = strfind(a0,[1 0]);
        numSigBins = endTrue - startTrue;
        use = find(numSigBins > 2);
        good = arrayfun(@(x) startTrue(use(x)):endTrue(use(x))-1, 1:length(use), ...
            'uniformoutput', false);
        sig_pre(cell2mat(good)) = true;


        %plot
        ax = nexttile;
        shadedErrorBar(x, clustPre_pref_Avg, clustPre_pref_Sem, 'lineprops', ...
            {'color', [0.65 0.65 0.65], 'linewidth', 1})
        hold on
        shadedErrorBar(x, clustPre_nonpref_Avg, clustPre_nonpref_Sem, 'lineprops', ...
            {'color', 'k', 'linewidth', 1})
        xlim([-0.5 1])
        % ylim([-0.5 1])
        yl = ylim;
        scatter(x(xvec(sig_pre)), yl(2), '*k')
        plot([0 0], [yl(1) yl(2)], '--k', 'linewidth', 0.25)
        hold on
        % set(gca, 'ytick', [])
        if c == nc
            set(gca, 'xtick', [-0.5, 0, 1])
            xlabel(labels{a})
        else
            set(gca, 'xtick', [])
        end
        set(gca, 'TickDir', 'out'); box off;
        ax.YRuler.TickLabelGapOffset = 1;
    

        %find 3 consecutive significant bins
        sig_post = false(1, length(xvec));
        a0 = [0 pExp_post_sig 0]; % add 0 at beginning and end to find pattern. also shifts index so start is the first true value
        startTrue = strfind(a0,[0 1]);
        endTrue = strfind(a0,[1 0]);
        numSigBins = endTrue - startTrue;
        use = find(numSigBins > 2);
        good = arrayfun(@(x) startTrue(use(x)):endTrue(use(x))-1, 1:length(use), ...
            'uniformoutput', false);
        sig_post(cell2mat(good)) = true;


        ax2 = nexttile;
        shadedErrorBar(x, clustPost_pref_Avg, clustPost_pref_Sem, 'lineprops', ...
            {'color', [0.65 0.65 0.65], 'linewidth', 1})
        hold on
        shadedErrorBar(x, clustPost_nonpref_Avg, clustPost_nonpref_Sem, 'lineprops', ...
            {'color', 'k', 'linewidth', 1})
        xlim([-0.5 1])
        ylim([yl])
        scatter(x(xvec(sig_post)), yl(2), '*k')
        plot([0 0], [yl(1) yl(2)], '--k', 'linewidth', 0.25)
        hold on
        set(gca, 'ytick', [])
        if c == nc
            set(gca, 'xtick', [-0.5, 0, 1])
            xlabel(labels{a})
        else
            set(gca, 'xtick', [])
        end
        set(gca, 'TickDir', 'out'); box off;
        ax2.YRuler.TickLabelGapOffset = 1;

    end
    
end
clim([-0.5 1.5])
colorbar(ax2, 'ticks', [-0.5, 0, 1.5])
set(gcf,'renderer','painter')
set(gcf, 'units', 'centimeters', 'position', [10 10 20 12])


