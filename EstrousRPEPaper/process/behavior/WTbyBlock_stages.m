function [WTbyBlock, pval_rats] =...
    WTbyBlock_stages(f_ratlist, RatBehaviorData)
%This function plots wait times by block for single rats and over the
%population. It then averages and normalizes the wait times by block to low the third reward per rat.
%Finally, it plots each rat's wait time by stage and block and evaluates
%the effect of stage on high/low wait times at the third reward with an
%ANOVA.

frats = length(f_ratlist);
%initialize variables
WTbyBlock = [];
Rewards = [1 2 3 4 5];
if grouping == 1
    cycle{1} = {'Proestrus', 'Estrus'};
    cycle{2} = {'Metestrus', 'Diestrus'};
    cyclenames = {'fertile','nonfertile'};
else
    cycle = {'Proestrus', 'Estrus', 'Metestrus', 'Diestrus'};
    cyclenames = cycle;
end
for rat = 1:frats

    ANOVA_waittimes = []; %populate variables for ANOVA
    ANOVA_blocktypes = []; %populate variables for ANOVA
    ANOVA_stages = []; %populate variables for ANOVA
        
    ratname = char(f_ratlist(rat));
    disp(['wt: ' ratname ' ' num2str(rat) ' out of ' num2str(frats)])

    stagedS = RatBehaviorData.(ratname);
    stagedS = stagedS.S;
    waittimeidx = cellfun(@(x) isfield(x, 'wait_time'), stagedS.pd);
    stagedS.pd = stagedS.pd(waittimeidx);

    for e = 1:length(cycle)

        stageidx = cellfun(@(x) logical(sum(strcmp(x.Stage, cycle{e}))), stagedS.pd); %find where it matches cycle stage(s)
        stageS.pd = stagedS.pd(stageidx);
        
        numsess = length(stageS.pd);

        %initialize variables to populate
        wt = cell(numsess, 1);
        bl = cell(numsess, 1);
        rewvol = cell(numsess, 1);
        high = cell(1, length(Rewards));
        low = cell(1, length(Rewards));
        mean_high = NaN(1, length(Rewards));
        sem_high = NaN(1, length(Rewards));
        mean_low = NaN(1, length(Rewards));
        sem_low = NaN(1, length(Rewards));

        for j = 1:numsess
            wts = stageS.pd{j}.wait_time;            
            bl{j, 1} = stageS.pd{j}.Block;
            rewvol{j, 1} = stageS.pd{j}.RewardAmount;
            opts = double(stageS.pd{j}.optout); %otherwise is logical
            delays = stageS.pd{j}.RewardDelay;

            %only keep wait times during opt outs
            wts(opts == 0 & delays ~= 100) = NaN;
            
            %remove multiple opt outs during catch trials in a row
            catch_index = find(delays == 100 & opts == 1);   
            d = [0; diff(catch_index)];
            wts(catch_index(d==1)) = NaN; 
            wt{j, 1} = wts;                   
           
        end

        waits = cell2mat(wt); %make one vector of wait times across sessions 
        waits(waits>(mean(waits, 'omitnan')+3*std(waits, 'omitnan'))) = NaN; %remove outliers
        blocktype = cell2mat(bl); %make one vector of blocks across sessions
        rewardvolume = cell2mat(rewvol); %make one vector of reward volume across sessions
        
        %z-score across all blocks
        waits_z = (waits-mean(waits, 'omitnan'))./std(waits, 'omitnan');

        %change reward volumes to 1, 2, 3, 4, or 5 (e.g. 1 = 4 or 5 ul, %etc.)
        rew_types = convertreward(rewardvolume);

        %find high block wait times (only catch trials)                 
        for rew = 1:length(Rewards)
            wtbyrew_high = waits(blocktype == 2 &...
                rew_types == Rewards(rew));
            high{rew} = wtbyrew_high;
            mean_high(1, rew) = mean(wtbyrew_high, 'omitnan');
            sem_high(1, rew) = std(wtbyrew_high, 'omitnan')./...
                sqrt(sum(~isnan(wtbyrew_high)));    
        end                  

        %find low block wait times (only catch trials)                       
        for rew = 1: length(Rewards)               
            wtbyrew_low = waits(blocktype == 3 &...
                rew_types == Rewards(rew));
            low{rew} = wtbyrew_low;
            mean_low(1, rew) = mean(wtbyrew_low, 'omitnan');
            sem_low(1, rew) = std(wtbyrew_low, 'omitnan')./...
                sqrt(sum(~isnan(wtbyrew_low))); 
        end                

        %divide averages by third reward during low block 
        mean_high_norm = mean_high/mean_low(3);
        mean_low_norm = mean_low/mean_low(3);
        
        if grouping==1
            cyclename = cell2mat(cycle{e});
        else
            cyclename = cycle{e};
        end

        delta = mean_low(3)-mean_high(3);
        
        WTbyBlock.(cyclenames{e}).mean_high_norm{rat} = mean_high_norm; %average wait time per reward for high block, normed to 3rd reward of low block
        WTbyBlock.(cyclenames{e}).mean_low_norm{rat} = mean_low_norm; %average wait time per reward for low block, normed to 3rd reward of low block            
        WTbyBlock.(cyclenames{e}).mean_high{rat} = mean_high; %average per reward for high block
        WTbyBlock.(cyclenames{e}).mean_low{rat} = mean_low; %average per reward for low block
        WTbyBlock.(cyclenames{e}).sem_high{rat} = sem_high; %SEM per reward for high block
        WTbyBlock.(cyclenames{e}).sem_low{rat} = sem_low; %SEM per reward for low block 
        WTbyBlock.(cyclenames{e}).waits{rat} = waits(rew_types == Rewards(3)); %wait times at 3rd reward
        WTbyBlock.(cyclenames{e}).waits_z{rat} = waits_z(rew_types == Rewards(3)); %wait times at 3rd reward
        WTbyBlock.(cyclenames{e}).block{rat} = blocktype(rew_types == Rewards(3)); %blocks at 3rd reward
        WTbyBlock.(cyclenames{e}).delta(rat) = delta;

        ANOVA_waittimes = [ANOVA_waittimes, waits(rew_types == Rewards(3))'];
        ANOVA_blocktypes = [ANOVA_blocktypes, blocktype(rew_types == Rewards(3))'];
        ANOVA_stages = [ANOVA_stages; repmat({cyclename}, length(waits(rew_types == Rewards(3))), 1)];
    
    end
    
    p_vals = anovan(ANOVA_waittimes(ANOVA_blocktypes~=1), {ANOVA_stages(ANOVA_blocktypes~=1)',...
        ANOVA_blocktypes(ANOVA_blocktypes~=1)},...
        'model', 'interaction', 'display','off');

    WTbyBlock.ANOVA_pvals{rat} = p_vals; %p-values for effects of block + stage + their interaction on wait time
    
end

%% STATS and PLOT
%%%plot all rats individually%%%
%%determine which rats have a significant effect of stage and which stage%%

if plot_arg==1

    if ismac
        figdir = '/Volumes/server/Estrus behavior figs/Wait Time/individual plots/grouped_stages/';
    elseif ispc
        figdir = 'Z:\Estrus behavior figs\Wait Time\individual plots\grouped_stages\';    
    end
    x = [1 2 3 4 5];
    cyclecolors = {[0.4940, 0.1840, 0.5560];[0 0 0];[0 1 1];[1 0 0]};
    wt_plots = zeros(2,1); % initialize plots so that they can have legends for block
    b_colors = [{[1 0 0]};{[0 0 1]}];
    
    pval_rats = NaN(frats, 1); %for ANOVA p-values by rat
    delta = NaN(frats, length(cycle)); %for mean high/mean low at 3rd reward
    
    for rat = 1:frats
        
        figure; 
        set(gcf, 'Color', [1 1 1], 'Units', 'Inches', 'Position', [0, 0, 17, 9],...
            'PaperUnits', 'Inches', 'PaperSize', [9, 5], 'renderer','painters')
    
        ratname = cell2mat(f_ratlist(rat));
    
        ratmax = NaN(2, 1);
        ratmin = NaN(2, 1);    
                
        for e = 1:length(cycle)
            subplot(1, length(cycle), e)
            
            mean_hi = WTbyBlock.(cyclenames{e}).mean_high{rat}; %mean high wts 
            mean_lo = WTbyBlock.(cyclenames{e}).mean_low{rat}; %mean low wts 
    
            er_hi = WTbyBlock.(cyclenames{e}).sem_high{rat}; %sem high wts
            er_lo = WTbyBlock.(cyclenames{e}).sem_low{rat}; %sem low wts
    
            %adaptation per stage
            delta(rat, e) = WTbyBlock.(cyclenames{e}).delta; %avg high@3rd/low@3rd reward
            
            %plot avg and sem
            if sum(~isnan([mean_hi mean_lo])) == 6
                wt_plots(1) = plot(x, mean_hi, 'Color',...
                    cell2mat(b_colors(1)), 'LineWidth', 2); hold on
                wt_plots(2) = plot(x, mean_lo, 'Color',...
                    cell2mat(b_colors(2)), 'LineWidth', 2); hold on    
    
                errorbar(x, mean_hi, er_hi,...
                    'Color', cell2mat(b_colors(1)), 'LineWidth', 2,...
                    'capsize', 0);
                errorbar(x, mean_lo, er_lo,...
                    'Color', cell2mat(b_colors(2)), 'LineWidth', 2,...
                    'capsize', 0);
    
                ratmax(e, 1) = max([mean_hi mean_lo]) + max([er_hi er_lo]);
                ratmin(e, 1) = min([mean_hi mean_lo]) - max([er_hi er_lo]);
    
                title(cyclenames{e})
                xticks([1 2 3 4 5])      
                xlim([0.5 5.5])
                ylabel('Wait Time (raw)')
                xlabel('Reward Level')
                grid off
                set(gca, 'TickDir', 'out'); box off
                axis square   
            end
        end
    
        if sum(~isnan([mean_hi mean_lo])) == 6
            for e = 1:length(cycle)
    
                ylimmax = max(ratmax);
                ylimmin = min(ratmin);
    
                subplot(1, length(cycle), e)
                ylim([ylimmin ylimmax])
    
            end
        end
    
        ANOVA_pvals = WTbyBlock.ANOVA_pvals{rat};
               
        pval_rats(rat) = ANOVA_pvals(3); %save interaction p-value for each rat
        
        sgtitle([ratname 'ANOVA'...
            'effect of stage: p =' num2str(round(ANOVA_pvals(1), 2))...
            'effect of block: p =' num2str(round(ANOVA_pvals(2), 2))...
            'interaction: p =' num2str(round(ANOVA_pvals(3), 2))]);
        
        saveas(gcf, strcat(figdir, ratname, '_ANOVA_', expt, '.png'))    
        saveas(gcf, strcat(figdir, ratname, '_ANOVA_', expt, '.pdf'))    
        saveas(gcf, strcat(figdir, ratname, '_ANOVA_', expt, '.fig'))  
            
        close all
    
    end

end

if frats > 1 && plot_arg==1

    %%%%%%%%%%% plot deltas %%%%%%%%%%%%%%%%%
    figdir = 'Z:\Estrus behavior figs\Wait Time\delta\individualpoints\';    
    figure; 
    set(gcf,'color','w','renderer','painters');
    plot([1 2], [mean(delta(:,1),'omitnan') mean(delta(:,2),'omitnan')],...
        '.k', MarkerSize=15); hold on
    errorbar([1 2], [mean(delta(:,1),'omitnan') mean(delta(:,2),'omitnan')],...
        [std(delta(:,1), 'omitnan')./sqrt(sum(~isnan(delta(:,1))))...
        std(delta(:,2), 'omitnan')./sqrt(sum(~isnan(delta(:,2))))],...
        color='k', linewidth=2, capsize=0); hold on
    xlim([0.5 2.5])
    for rat=1:length(f_ratlist)
        plot([1 2], [delta(rat,1) delta(rat,2)], '-',...
            color=[0 0 0 0.25]); hold on
    end
    pval = signrank(delta(:,1), delta(:,2));
    subtitle(['p=' num2str(pval)])
    xticks(1:2)
    xticklabels({'fertile','non-fertile'})
    ylabel('Low-high at 16 ul')
    axis square
    grid off
    set(gca, 'TickDir', 'out'); box off
    saveas(gcf, strcat(figdir, 'wtbyblock_delta_', expt, '.png'))
    saveas(gcf, strcat(figdir, 'wtbyblock_delta_', expt, '.pdf'))
    saveas(gcf, strcat(figdir, 'wtbyblock_delta_', expt, '.fig'))

    %%%%%%%%%%% plot wait times by block and stage %%%%%%%%%%%
    wt_plots = zeros(2,1); % initialize plots so that they can have legends for block
    b_colors = [{[1 0 0]};{[0 0 1]}];
    x = [1 2 3 4 5];

    figure; 
    set(gcf,'color','w','renderer','painters');

    for e = 1:length(cycle)

        hi = NaN(length(frats), length(x));
        lo = NaN(length(frats), length(x));

        subplot(1, length(cycle), e)

        %plot individual rats
        for rat = 1: frats     

%             hi(rat, :) = WTbyBlock.(cyclenames{e}).mean_high_norm{rat}; %mean high wts per rat
            hi(rat, :) = WTbyBlock.(cyclenames{e}).mean_high{rat}; %mean high wts per rat
%             lo(rat, :) = WTbyBlock.(cyclenames{e}).mean_low_norm{rat}; %mean low wts per rat
            lo(rat, :) = WTbyBlock.(cyclenames{e}).mean_low{rat}; %mean low wts per rat

            plot(x, hi(rat, :), 'Color', [cell2mat(b_colors(1)), 0.3],...
                'LineWidth', 1); hold on
            plot(x, lo(rat, :), 'Color', [cell2mat(b_colors(2)), 0.3],...
                'LineWidth', 1); hold on

        end

        %plot avg
        wt_plots(1, 1) = plot(x, mean(hi, 'omitnan'), 'Color', cell2mat(b_colors(1)),...
            'LineWidth', 2); hold on
        wt_plots(2, 1) = plot(x, mean(lo, 'omitnan'), 'Color', cell2mat(b_colors(2)),...
            'LineWidth', 2); hold on    

%         ylim([0.4 1.5])
        title(cyclenames{e})
        ylabel('Wait Time High/Low Block @3rd Reward ')
        xlabel('Reward Level')
        grid off
        set(gca, 'TickDir', 'out'); box off
        axis square    
    end

    %%%plot average and sem%%%
    figure; 
    set(gcf,'color','w','renderer','painters');

    if ismac
        figdir = '/Volumes/server/Estrus behavior figs/Wait Time/by block/summary/';
    elseif ispc
        figdir = 'Z:\Estrus behavior figs\Wait Time\by block\summary\';    
    end

    ratsmax = NaN(2, 1);
    ratsmin = NaN(2, 1);   

    for e = 1:length(cycle)

        high_rats = NaN(1, length(frats));
        low_rats = NaN(1, length(frats));

        subplot(1, length(cycle), e)
        for rat = 1: frats     
            wts = WTbyBlock.(cyclenames{e}).waits_z{rat}; %wait times at 3rd reward
            bls = WTbyBlock.(cyclenames{e}).block{rat}; %blocks at 3rd reward
            high_rats(1, rat) = mean(wts(bls==2), 'omitnan');
            low_rats(1, rat) = mean(wts(bls==3), 'omitnan');
        end

        %plot avg and sem
        meanhi = mean(high_rats, 'omitnan');
        meanlo = mean(low_rats, 'omitnan');
        er_hi = std(high_rats, 'omitnan');
        er_lo = std(low_rats, 'omitnan');

        plot(1, meanhi, 'Color', cell2mat(b_colors(1)),...
            'LineWidth', 2); hold on
        plot(2, meanlo, 'Color', cell2mat(b_colors(2)),...
            'LineWidth', 2); hold on    

        errorbar(1, meanhi, er_hi,...
            'Color', cell2mat(b_colors(1)), 'LineWidth', 2, 'capsize', 0); hold on   
        errorbar(2, meanlo, er_lo,...
            'Color', cell2mat(b_colors(2)), 'LineWidth', 2, 'capsize', 0); hold on   

        %find number of rats included in figure
        N = sum(~isnan(low_rats));

        title(strcat(cyclenames{e}, ':', 'N = ', num2str(N)))
        xticks([1 2])        
        xlim([0.5 2.5])
        ylabel('Wait Time at 16 ul (z-scored, s)')
        xticklabels({'High', 'Low'})
        grid off
        set(gca, 'TickDir', 'out'); box off
        axis square    

        ratsmax(e, 1) = max([meanhi meanlo]) + max([er_hi er_lo]);
        ratsmin(e, 1) = min([meanhi meanlo]) - max([er_hi er_lo]);

    end 

    ylimmax = max(ratsmax);
    ylimmin = min(ratsmin);

    for e = 1:length(cycle)
        subplot(1, length(cycle), e)
        ylim([ylimmin ylimmax])
    end

    saveas(gcf, strcat(figdir, 'wtbyblock_at16ul_summary_minsess', num2str(minsess),...
        '_', expt, '.png'))
    saveas(gcf, strcat(figdir, 'wtbyblock_at16ul_summary_minsess', num2str(minsess),...
        '_', expt, '.pdf'))
    saveas(gcf, strcat(figdir, 'wtbyblock_at16ul_summary_minsess', num2str(minsess),...
        '_', expt, '.fig'))
    
    %%%plot average and sem normed%%%
    wt_plots = zeros(2,1); % initialize plots so that they can have legends for block

    figure; 
    set(gcf,'color','w','renderer','painters');

    if ismac
        figdir = '/Volumes/server/Estrus behavior figs/Wait Time/by block/summary/';
    elseif ispc
        figdir = 'Z:\Estrus behavior figs\Wait Time\by block\summary\';    
    end

    for e = 1:length(cycle)

        hi = NaN(length(frats), length(x));
        lo = NaN(length(frats), length(x));

        linetypes = {'-', '--'};

        for rat = 1: frats     
%             hi(rat, :) = WTbyBlock.(cyclenames{e}).mean_high_norm{rat}; %mean high wts per rat
            hi(rat, :) = WTbyBlock.(cyclenames{e}).mean_high{rat}; %mean high wts per rat
%             lo(rat, :) = WTbyBlock.(cyclenames{e}).mean_low_norm{rat}; %mean low wts per rat
            lo(rat, :) = WTbyBlock.(cyclenames{e}).mean_low{rat}; %mean low wts per rat
        end

        %plot avg and sem
        wt_plots(1) = plot(x, mean(hi, 'omitnan'), linetypes{e},...
            'Color', cell2mat(b_colors(1)),...
            'LineWidth', 2); hold on
        wt_plots(2) = plot(x, mean(lo, 'omitnan'), linetypes{e},...
            'Color', cell2mat(b_colors(2)),...
            'LineWidth', 2); hold on    

        errorbar(x, mean(hi, 'omitnan'), std(hi, 'omitnan')./sqrt(sum(~isnan(hi))),...
            linetypes{e}, 'Color', cell2mat(b_colors(1)), 'LineWidth', 2, 'capsize', 0); hold on   
        errorbar(x, mean(lo, 'omitnan'), std(lo, 'omitnan')./sqrt(sum(~isnan(lo))),...
            linetypes{e}, 'Color', cell2mat(b_colors(2)), 'LineWidth', 2, 'capsize', 0); hold on   

        %find number of rats included in figure
        N = sum(~isnan(lo(:, 3)));

    %     legend(wt_plots, 'High Block', 'Low Block')
        xticks([1 2 3 4 5])        
        xlim([0.5 5.5])
%         ylim([0.85 1.1])
%         ylim([10 13])
%         ylabel('High/Low Block @3rd Reward ')
        ylabel('Wait time (s)')
        xlabel('Reward Level')
        grid off
        set(gca, 'TickDir', 'out'); box off
        axis square    
    end 
    title(['N=' num2str(N)])

    saveas(gcf, strcat(figdir, 'wtbyblock_summary_minsess', num2str(minsess),...
        '_', expt, '.png'))
    saveas(gcf, strcat(figdir, 'wtbyblock_summary_minsess', num2str(minsess),...
        '_', expt, '.pdf'))
    saveas(gcf, strcat(figdir, 'wtbyblock_summary_minsess', num2str(minsess),...
        '_', expt, '.fig'))

    %%%plot average and sem normed%%%
    delta_plots = zeros(2,1); % initialize plots so that they can have legends for block

    figure; 
    set(gcf,'color','w','renderer','painters');

    if ismac
        figdir = '/Volumes/server/Estrus behavior figs/Wait Time/by block/deltas/';
    elseif ispc
        figdir = 'Z:\Estrus behavior figs\Wait Time\by block\deltas\';    
    end

    deltas = NaN(length(frats), 2);

    for e = 1:length(cycle)

        for rat = 1:frats     
            hi = WTbyBlock.(cyclenames{e}).mean_high{rat}; %mean high wts per rat
            lo = WTbyBlock.(cyclenames{e}).mean_low{rat}; %mean low wts per rat
            deltas(rat, e) = lo(:, 3) - hi(:, 3);
        end
        %plot avg and sem
        delta_plots(e, 1) = plot(e, mean(deltas(:, e), 'omitnan'),...
            'Color', cell2mat(b_colors(1)),...
            'LineWidth', 2); hold on 
        errorbar(e, mean(deltas(:, e), 'omitnan'),...
            std(deltas(:, e), 'omitnan')./sqrt(sum(~isnan(deltas(:, e)))),...
            'Color', cyclecolors{e}, 'LineWidth', 2, 'capsize', 0); hold on   
    end 
    for rat = 1:frats     
        %plot individual lines
        plot(1:length(cycle), deltas(rat, :),...
            'Color', [0 0 0 0.3], 'LineWidth', 1); hold on 
    end
    xticks(1:length(cycle))   
    if grouping==1
        thesecyclenames{1} = cell2mat(cycle{1});
        thesecyclenames{2} = cell2mat(cycle{2});
        pval = signrank(deltas(:, 1), deltas(:, 2));
        subtitle(['signrank p=' num2str(pval)])
    else
        thesecyclenames = cycle;
        pval = kruskalwallis(deltas(:),...
            [ones(frats,1); 2*ones(frats,1);...
            3*ones(frats,1); 4*ones(frats,1)], 'off');
        subtitle(['kruskal wallis p=' num2str(pval)])
    end
    xticklabels(thesecyclenames)
    xlim([0.5 length(cycle)+0.5])
    ylabel('Wait time at 16 ul low - high (s)')
    grid off
    set(gca, 'TickDir', 'out'); box off
    axis square    
    
    N = sum(~isnan(deltas(:, 1)));
    title(strcat('N=', num2str(N)))
    

    saveas(gcf, strcat(figdir, 'wtbyblock_deltaindiv_summary_minsess', num2str(minsess),...
        '_', expt, '.png'))
    saveas(gcf, strcat(figdir, 'wtbyblock_deltaindiv_summary_minsess', num2str(minsess),...
        '_', expt, '.pdf'))
    saveas(gcf, strcat(figdir, 'wtbyblock_deltaindiv_summary_minsess', num2str(minsess),...
        '_', expt, '.fig'))

    %%%plot deltas as individual points and with summary statistics%%%
    delta_plots = zeros(length(cycle),1); % initialize plots so that they can have legends for estradiol level
    x = 1:length(cycle);

    figure; 
    set(gcf,'color','w','renderer','painters');
    if ismac
%         figdir = '/Volumes/server/Estrus behavior figs/Wait Time/delta/individualpoints/';
    elseif ispc
        figdir = 'Z:\Estrus behavior figs\Wait Time\delta\individualpoints\';    
    end
    
    %plot individual points
    for rat = 1: frats   
        ratdeltas = [];
        for e = 1:length(cycle)
            jitx = randn(1)*0.05;
            plot(e+jitx, delta(rat, e), 'o', 'Color',...
                [cyclecolors{e}, 0.3],...
                'LineWidth', 1); hold on %block individual data points
            ratdeltas = [ratdeltas, delta(rat, e)];
        end
        
        plot(x, ratdeltas, '-', 'Color', [0.5 0.5 0.5],...
            'LineWidth', 1); hold on %block connecting line between rats
    end

    %plot average:
    %connecting line
    mean_delta_rats = NaN(length(cycle), 1);
    for e = 1:length(cycle)
        mean_delta_rats(e) = mean(delta(:, e), 'omitnan');
    end
    plot(x, mean_delta_rats, '-k', 'LineWidth', 2); hold on   
    for e = 1:length(cycle)
        delta_plots(e, 1) = plot(e, mean(delta(:, e), 'omitnan'), '+',...
            'Color', cyclecolors{e}, 'LineWidth', 2); hold on   
    end

    %plot average and sem 
    for e = 1:length(cycle)
        delta_plots(e, 1) = plot(e, mean(delta(:, e), 'omitnan'), 'o',...
            'Color', cyclecolors{e}, 'LineWidth', 2); hold on
        errorbar(e, mean(delta(:, e), 'omitnan'),...
            std(delta(:, e), 'omitnan')./sqrt(sum(~isnan(delta(:, e)))),...
            'Color', cyclecolors{e}, 'LineWidth', 2, 'capsize', 0); hold on        
    end

    %delta_p = friedman(delta,11);

    xlim([0 length(cycle)+1]);
    xticks(x)
    if grouping==1
        cyclenames{1} = cell2mat(cycle{1});
        cyclenames{2} = cell2mat(cycle{2});
    else
        cyclenames = cycle;
    end
    xticklabels(cyclenames)
    ylabel('Wait Time @3rd Reward (High/Low Blocks)')
    grid off
    set(gca, 'TickDir', 'out'); box off
    axis square

    saveas(gcf, strcat(figdir, 'minsess_', num2str(minsess),...
        '_', expt, '.png'))
    saveas(gcf, strcat(figdir, 'minsess_', num2str(minsess),...
        '_', expt, '.pdf'))

    %%%plot delta summary statistics%%%
    delta_plots = zeros(length(cycle),1); % initialize plots so that they can have legends for estradiol level
    x = 1:length(cycle);

    figure; 
    set(gcf,'color','w','renderer','painters');
    if ismac
        figdir = '/Volumes/server/Estrus behavior figs/Wait Time/delta/summary/';
    elseif ispc
        figdir = 'Z:\Estrus behavior figs\Wait Time\delta\summary\';    
    end
    %plot average and sem 
    for e = 1:length(cycle)
        delta_plots(e, 1) = plot(e, mean(delta(:, e), 'omitnan'), 'o',...
            color=cyclecolors{e}, LineWidth=1); hold on
        errorbar(e, mean(delta(:, e), 'omitnan'),...
            std(delta(:, e), 'omitnan')./sqrt(sum(~isnan(delta(:, e)))),...
            Color=cyclecolors{e}, LineWidth=1, capsize=0); hold on        
    end

    if grouping == 1
        delta_p = signrank(delta(:, 1), delta(:, 2));
    else
        delta_p = kruskalwallis(delta);
    end
    xlim([0 length(cycle)+1]);
    xticks(x)
    xticklabels(cyclenames)
    ylabel('Wait Time @3rd Reward (1 - High/Low Blocks)')
    title({strcat('Estrous stage:'),...
        strcat('p = ', " ", num2str(delta_p))});
%     ymin = min(mean(delta(:), 'omitnan'))-std(delta(:), 'omitnan');
%     ylim([ymin 1])
%     ymax = max(mean(delta(:), 'omitnan'))+std(delta(:), 'omitnan');
%     ylim([ymin ymax])
    grid off
    set(gca, 'TickDir', 'out'); box off
    axis square

    saveas(gcf, strcat(figdir, 'wtbyblock_delta_minsess_', num2str(minsess),...
        '_', expt, '.png'))
    saveas(gcf, strcat(figdir, 'wtbyblock_delta_minsess_', num2str(minsess),...
        '_', expt, '.pdf'))
    saveas(gcf, strcat(figdir, 'wtbyblock_delta_minsess_', num2str(minsess),...
        '_', expt, '.fig'))

    %%%%%%%%%%% plot histogram of p-values %%%%%%%%%%%%%%%%%
    if ismac
        figdir = '/Volumes/server/Estrus behavior figs/Wait Time/pval_stageeffect/';
    elseif ispc
        figdir = 'Z:\Estrus behavior figs\Wait Time\pval_stageeffect\';    
    end    
    figure;
    histogram(pval_rats, binwidth=0.01, facecolor='k'); hold on
    xline([0.05 0.05], '--k')
    grid off
    set(gca, 'TickDir', 'out'); box off
    axis square   
    title([num2str(sum(pval_rats<0.05)) ' significant rats'])
    saveas(gcf, [figdir, 'pvalhist_', expt, '.png'])
    saveas(gcf, [figdir, 'pvalhist_', expt, '.pdf'])

end

end
