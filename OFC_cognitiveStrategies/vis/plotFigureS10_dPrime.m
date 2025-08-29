function plotFigureS10_dPrime(processedEphysPath)
% Plot Figure S10. Plot fraction of block significant cells that go into 
% incongruent response plots for each recording area. Plot discriminability 
% index for large vs small volumes for each recording area

% YOU MUST RUN processEphysData AND processDPrime FOR OFC,
% MOTOR CORTEX, AND PIRIFORM CORTEX DATA FOR BOTH EXPERT AND NAIVE RATS
% BEFORE RUNNING THIS FUNCTION. 

% INPUTS
%   processedEphysPath = local path to processed incongruent response and
%       d' data. For best results, processed data should be saved into folders 
%       structured as e.g. processedEphysPath/Expert/OFC. See below for
%       expected structure and naming. These can be changed if desired


%% Paths - You may need to change these to match how you have saved the 
% processed incongruent response and d' data for each area. Function
% assumes you have saved incongruent response and d' data for the same
% group and area (e.g. expert OFC) to the same folder.

OFCpath_exp = [processedEphysPath filesep 'Expert' filesep 'OFC'];
MotorPath_exp = [processedEphysPath filesep 'Expert' filesep 'Motor'];
PiriformPath_exp = [processedEphysPath filesep 'Expert' filesep 'Piriform'];

expertPaths = {OFCpath_exp MotorPath_exp PiriformPath_exp};


OFCpath_naive = [processedEphysPath filesep 'Naive' filesep 'OFC'];
MotorPath_naive = [processedEphysPath filesep 'Naive' filesep 'Motor'];
PiriformPath_naive = [processedEphysPath filesep 'Naive' filesep 'Piriform'];

naivePaths = {OFCpath_naive MotorPath_naive PiriformPath_naive};


%% general functions
sem = @(x) std(x,'omitnan')./sqrt(size(x,1));
setDefaultFigProps
fsize = [10 8 20 14];

alignto = {'COFF' 'SON' 'SOFF' 'Rew'};
ne = length(alignto);

wndw = [-0.5 4];
areas = {'OFC' 'Motor cortex' 'Piriform cortex'};

%%
figure; hold on
tiledlayout(length(expertPaths), ne+1, 'tilespacing', 'compact')
    
for a = 1:length(expertPaths)

    incongruent_exp = load([expertPaths{a} filesep, 'incongruentResponses.mat']);

    incongruent_naive = load([naivePaths{a}, filesep, 'incongruentResponses.mat']);

    d_exp = load([expertPaths{a}, filesep, 'd.mat']);

    d_naive = load([naivePaths{a}, filesep, 'd.mat']);

    
    % Get fraction of block significant cells for each alignment
    sig_exp = arrayfun(@(x) sum(incongruent_exp.info.(alignto{x})), 1:ne);
    fracSig_exp = sig_exp ./ height(incongruent_exp.info);
    
    sig_naive = arrayfun(@(x) sum(incongruent_naive.info.(alignto{x})), 1:ne);
    fracSig_naive = sig_naive ./ height(incongruent_naive.info);
    
    
    for jj = 1:ne
        %binomial confidence intervals
        [~, pci_exp(jj,:)] = binofit(sig_exp(jj), height(incongruent_exp.info));
        [~, pci_naive(jj,:)] = binofit(sig_naive(jj), height(incongruent_naive.info));
    
        %remove any infs from d'
        d_exp.d_vol{jj}(isinf(d_exp.d_vol{jj})) = nan;
        d_naive.d_vol{jj}(isinf(d_naive.d_vol{jj})) = nan;
    end

    %Average d'
    meansE_vol = arrayfun(@(x) mean(d_exp.d_vol{x}, 'omitnan'), ...
        1:length(alignto), 'uniformoutput', false);
    semE_vol = arrayfun(@(x) sem(d_exp.d_vol{x}), 1:length(alignto), ...
        'uniformoutput', false);

    meansN_vol = arrayfun(@(x) mean(d_naive.d_vol{x}, 'omitnan'), ...
        1:length(alignto), 'uniformoutput', false);
    semN_vol = arrayfun(@(x) sem(d_naive.d_vol{x}), 1:length(alignto), ...
        'uniformoutput', false);


    % plot fraction of significant cells
    x = 1:4;
    y = [fracSig_exp; fracSig_naive];
    err_exp = [fracSig_exp - pci_exp(:,1)'; pci_exp(:,2)' - fracSig_exp];
    err_naive = [fracSig_naive - pci_naive(:,1)'; pci_naive(:,2)' - fracSig_naive];
    err_pos = [err_exp(1,:); err_naive(1,:)];
    err_neg = [err_exp(2,:); err_naive(2,:)];

    [ngroups, nbars] = size(y);
   

    nexttile
    b = bar(x, y);
    hold on
    for ii = 1:ngroups
        x2 = b(ii).XEndPoints; %midpoint of each bar for group
        errorbar(x2, y(ii,:), err_pos(ii,:), err_neg(ii,:), '.k');
    end
    b(1).FaceColor = 'k'; b(1).FaceAlpha = 0.6;
    b(2).FaceColor = 'm'; b(2).FaceAlpha = 0.6;
    ylim([0 0.5])
    xticks(1:4)

    yticks(0:0.1:0.5)
    ylabel('Fraction')
    set(gca, 'TickDir', 'out'); box off;
    % axis square
    ax1 = gca;
    ax1.YRuler.TickLabelGapOffset = 1;
    if a == 1
        legend({'Expert' 'Naive'})
    end
    if a == length(expertPaths)
        xticklabels(alignto)
    else
        set(gca, 'XTickLabel', []);
    end
    title(areas{a})


    for ii = 1:ne
        xvec = wndw(1):0.05:wndw(2);

        ax(ii) = nexttile;
        shadedErrorBar(xvec, meansE_vol{ii}, semE_vol{ii})
        shadedErrorBar(xvec, meansN_vol{ii}, semN_vol{ii}, 'lineprops', 'm')
        hold on
        xline(0, '--k', 'LineWidth', 1)
        ylim([0 0.45])
        xlim([-0.5 2])
        set(gca, 'TickDir', 'out'); box off
        % axis square
        if ii ==1
            ylabel("volume d'")
        else
            set(gca,'YTickLabel',[]);
        end
        if a ~= length(expertPaths)
            set(gca, 'XTickLabel', []);
        end
    end

end
set(gcf, 'Renderer', 'painters')
set(gcf, 'units', 'centimeters', 'position', fsize)