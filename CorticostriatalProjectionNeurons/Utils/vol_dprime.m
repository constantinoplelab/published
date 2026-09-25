function [dprime, sig,goods,dvar,dprime_sig_cells] = vol_dprime(trials,goods,all_S,all_SU,all_index,event,win,varargin)
% compares the d'prime for individual cells, and mean, compares the
% mean to the shuffle data,and shows the mean dprime with the shuffle
% subtracted

%inputs: 
    % trials - cells x 2 for each condition
    % goods - cells to consider
    % all_S - S structs
    % all_SU - SU structs
    % all_index - index of all cells
    % event - event you want to consider
    % win - time window to display
    % varargin - 
    %   1. "signed"/"unsigned"
    %   2. # of interations to average of the shuffle
    %   3. plot or not
 %output: 
    % dprime of the data - if you provide the number of interations for the
    %   shuffle, dprime is the shuffle subtracted dprime
    % sig - time points that the mean dprime is larger than the shuffle 95th percentile
    % dvar - variance that is used in d' calculation
    % dprime_sig_cells = d prime for individual clusters with values
    %   greater than the 95th percentile

if length(varargin)>2
    PlotIt = varargin{3};
end

% remove sessions with <4 trials in either category. 
trialnum = cellfun(@length,trials);
trialnum = min(trialnum,[],2);
trials(trialnum<4,:) = [];
if sum(trialnum<4)>0
goods(trialnum<4) =[];
end

xvec = all_SU{goods(1)}.xvec.COFF;

% d' for data
if strcmp(varargin{1},'signed')
    [d,dvar,dprime,~] = makeDprime(trials, goods, all_S,all_SU,all_index,event,0);
else 
    [d,dvar,dprime,~] = makeDprime(trials, goods, all_S,all_SU,all_index,event,0,1);
end
% remove inf from the dprime
% inf arrises from only have 1 trial, so variance is impossible to
% calculate
dprime(abs(dprime)>1000) = NaN;

% shuffle d' for comparison
for i =1:numel(trials)
    trials{i} = trials{i}';
end

if nargin>8
    %% single and multiple shuffle
    % % single shuffle
    % shuffle = makeshuffle_balanced(trials);
    % shuffle = shuffle.COFF;
    % if strcmp(varargin{1},'signed')
    % [d_shuffle_s,dvar_shuffle,dprime_shuffle_s,h_shuffle] = makeDprime(shuffle, goods, all_S,all_SU,all_index,event);
    % else 
    % [d_shuffle_s,dvar_shuffle,dprime_shuffle_s,h_shuffle] = makeDprime(shuffle, goods, all_S,all_SU,all_index,event,1);
    % end
    % dprime_shuffle_s(abs(dprime_shuffle_s)>1000) = NaN;

    % multiple shuffle interations
    for s = 1:varargin{2}
    shuffle = makeshuffle_balanced(trials,s);
    shuffle = shuffle.COFF;
    if strcmp(varargin{1},'signed')
    [d_shuffle(:,:,s),dvar_shuffle,dprime_shuffle(:,:,s),h_shuffle] = makeDprime(shuffle, goods, all_S,all_SU,all_index,event,0);
    else 
    [d_shuffle(:,:,s),dvar_shuffle,dprime_shuffle(:,:,s),h_shuffle] = makeDprime(shuffle, goods, all_S,all_SU,all_index,event,0,1);
    end
    end %interations of shuffle
    
    % remove Inf before any averaging
    dprime_shuffle(abs(dprime_shuffle)>1000) = NaN;
    
    %average
    dprime_shuffle_avg = mean(dprime_shuffle,3,'omitnan');  %across shuffle iterations
    dprime_shuffle_avg_shuffles = squeeze(mean(dprime_shuffle,1,'omitnan'))'; %across cell iterations
    
    % %plot the shuffle
    % figure
    % plot(xvec,dprime_shuffle_avg)

    % shuffle subtract the real dprime
    dprimeShufRemove = dprime-dprime_shuffle_avg;

    % Modify shuffle subtracted mean to only include significant values
    % Significant = greater than the 95 percentile of the shuffle mean
    % didn't use std bc we don't know that the shuffle will be normally
    % distributed.
    % this also currently only does greater - that could be a problem for
    % signed d'
    % also if you comment this part out, make sure you change the dprime
    % output
    dprime_shuffle_avg_quant = quantile(dprime_shuffle_avg,0.95, 1);
    dprime_shuffle_avg_shuffles_quant = quantile(dprime_shuffle_avg_shuffles,0.95, 1);


        %for the individual cell remove values less than the 95th
        %percentile
        dprime_ShuffRemove_cell = dprime-dprime_shuffle_avg;

    

    %make nan's zero for subtraction that doesn't result in everything NaN
    % dprime_shuffle_avg_zero=dprime_shuffle_avg;
    % dprime_shuffle_avg_zero(isnan(dprime_shuffle_avg))=0;
    % dprime_zero=dprime;
    % dprime_zero(isnan(dprime)) = 0;
    % 
    % dprimeShufRemove = dprime_zero-dprime_shuffle_avg_zero;
    

    if nargin<10 %plot or not
        dprime_data_avg = mean(dprime,'omitnan');
            dprime_shuffle_avg_avg = mean(dprime_shuffle_avg,'omitnan');

        sig = dprime_data_avg>dprime_shuffle_avg_shuffles_quant  |  dprime_data_avg<dprime_shuffle_avg_avg-dprime_shuffle_avg_shuffles_quant ;

        dprime = dprimeShufRemove;
    else


    figure('color','white')
    T = find(xvec>=win(1),1):find(xvec<=win(2),1,'last');
    % Fig1. Grey = each cell d'. Raw d'. 
    %       Bold = mean d'
    subplot(1,4,1)
    plot(xvec(T),dprime(:,T),'color',[0.5, 0.5, 0.5]); hold on;
    plot(xvec(T),mean(dprime(:,T),'omitnan')','k','LineWidth',2)
    box off
    set(gca,'TickDir','out')
    title('individual cell dprime ')
    xlabel('time (s)')
    ylabel('dprime')

    % Fig2. Grey = mean+sem of d' shuffle. d' shuffle = average all shuffle
    % iterations of each cell to get 1 d' psth for each cell. Then average
    % those psths to get population mean d' psth. 
    % Blue = mean+sem of d'.
    subplot(1,4,2)
    shadedErrorBar(xvec(T),mean(dprime(:,T),'omitnan'),sem(dprime(:,T),'omitnan'),'LineProps',{'color',[0 0 1]})
    shadedErrorBar(xvec(T),mean(dprime_shuffle_avg(:,T),'omitnan'),sem(dprime_shuffle_avg(:,T),'omitnan'),'LineProps',{'color',[0.5 0.5 0.5]})
    %shadedErrorBar(xvec(T),mean(dprime_shuffle_avg_shuffles(:,T),'omitnan'),sem(dprime_shuffle_avg_shuffles(:,T),'omitnan'),'LineProps',{'color',[1 0.5 0.5]})
    xline(0,'--')
    [t_test.h,t_test.p] = ttest(dprime(:,T)-dprime_shuffle_avg(:,T));
    %[h,p,ci,stats] = ttest2(dprime(:,T),dprime_shuffle_avg(:,T)); %should
    %be the same as the ttest above
    hold on;
    plot(xvec(T(t_test.h==1)),t_test.h(t_test.h==1),'*k')
    xlabel('time (s)')
    ylabel('dprime')
    legend({'data','shuffle'})
    title('Mean dprime')

    % Fig3. Grey = mean d' shuffle+95% confidence intervals. d' shuffle =
    % average all cells in one iteration to get 1 population d' psth for
    % each iteration. Then average those psths to get a population mean d'
    % psth. 
    % Blue = mean+sem d' data. 
    subplot(1,4,3)
    hold on;
    % get the average of the shuffle distribution
    dprime_shuffle_avg_avg = mean(dprime_shuffle_avg,'omitnan');
    xline(0,'--')
    
    %% use the distribution confidence intervals
    % %percentile takes the mean into account, so you need to subtract it
    % shadedErrorBar(xvec(T), dprime_shuffle_avg_avg,  dprime_shuffle_avg_quant(:,T)-dprime_shuffle_avg_avg, 'lineprops',{'color',[0 0 0]})
    % %plot(xvec(T),dprime_shuffle_avg_quant(:,T),'color',[0.5 0.5 0.5]);
    shadedErrorBar(xvec(T), dprime_shuffle_avg_avg(:,T),  dprime_shuffle_avg_shuffles_quant(:,T)-dprime_shuffle_avg_avg(:,T), 'lineprops',{'color',[0 0 0]})
    %plot(xvec(T),dprime_shuffle_avg_shuffles_quant(:,T),'color',[1 0.5 0.5]);
    
    % %% uses the distribution std
    % shadedErrorBar(xvec(T),dprime_shuffle_avg_avg,2*std(dprime_shuffle_avg(:,T),'omitnan'),'lineprops',{'color',[0 0.5 0]});
    % shadedErrorBar(xvec(T),dprime_shuffle_avg_avg,2*std(dprime_shuffle_avg_shuffles(:,T),'omitnan'),'lineprops',{'color',[0 0.5 0]});

    %data
    dprime_data_avg = mean(dprime,'omitnan');
    plot(xvec(T),dprime_data_avg(:,T),'color',[0 0 1])
    %legend({'across shuffle confidence','data','sig data'})
    title('compare to shuffle confidence intervals')
    xlabel('Time (s)')
    %bold regions where the data mean>confidence inter
    
    sig = dprime_data_avg>dprime_shuffle_avg_shuffles_quant  |  dprime_data_avg<dprime_shuffle_avg_avg-dprime_shuffle_avg_shuffles_quant ;
    dprime_data_avg2 = dprime_data_avg;
    dprime_data_avg2(sig==0) = NaN;
    plot(xvec(T),dprime_data_avg2(:,T),'color',[0 0 1],'LineWidth',2)

    % Fig4. Blue: Shuffle subtracted d' mean+sem. shuffle subtracted =raw
    % d' - shuffle averaged across iterations
    subplot(1,4,4)
    hold on;
    dprimeShufRemove2 = dprimeShufRemove;
    dprimeShufRemove2(:,sig==0) = NaN;
    shadedErrorBar(xvec(T),mean(dprimeShufRemove(:,T),'omitnan'),sem(dprimeShufRemove(:,T),'omitnan'),'LineProps',{'color',[0.1 0.1 1],'linewidth',1})
    plot(xvec(T),mean(dprimeShufRemove2(:,T),'omitnan'),'color',[0.1 0.1 1],'linewidth',2)

    % also include the significant values only for the shuffle-subtracted
    %shadedErrorBar(xvec(T),mean(dprime_ShuffRemove_SigOnly(:,T),'omitnan'),sem(dprime_ShuffRemove_SigOnly(:,T),'omitnan'),'LineProps',{'color',[0 0 0]});
    
    xline(0,'--')
    xlabel('time (s)')
    ylabel('dprime')
    title('Mean dprime shuffle subtracted')
    sgtitle(['dprime at ', event])
    end

dprime = dprimeShufRemove;
%dprime = dprime_ShuffRemove_SigOnly;

else 
    %% single shuffle only
    % single shuffle
    shuffle = makeshuffle_balanced(trials);
    shuffle = shuffle.COFF;
    if strcmp(varargin{1},'signed')
    [d_shuffle_s,dvar_shuffle,dprime_shuffle_s,h_shuffle] = makeDprime(shuffle, goods, all_S,all_SU,all_index,event);
    else 
    [d_shuffle_s,dvar_shuffle,dprime_shuffle_s,h_shuffle] = makeDprime(shuffle, goods, all_S,all_SU,all_index,event,1);
    end
    dprime_shuffle_s(dprime_shuffle_s>1000) = NaN;


    if nargin>9
    else
    %plot 
    figure
    T = find(xvec>=win(1),1):find(xvec<=win(2),1,'last');

    subplot(1,2,1)
    plot(xvec(T),dprime(:,T),'color',[0.5, 0.5, 0.5]); hold on;
    plot(xvec(T),mean(dprime(:,T),'omitnan')','k','LineWidth',2)
    box off
    set(gca,'TickDir','out')
    title('individual cell dprime ')
    xlabel('time (s)')
    ylabel('dprime')

    subplot(1,2,2)
    shadedErrorBar(xvec(T),mean(dprime(:,T),'omitnan'),sem(dprime(:,T),'omitnan'),'LineProps',{'color',[0 0 0]})
    shadedErrorBar(xvec(T),mean(dprime_shuffle_s(:,T),'omitnan'),sem(dprime_shuffle_s(:,T),'omitnan'),'LineProps',{'color',[0.5 0.5 0.5]})
    xline(0,'--')
    xlabel('time (s)')
    ylabel('mean dprime')
    legend({'data','shuffle'})
    title('Mean dprime')
    set(gca,'TickDir','out')

    end
    sig = NaN;

end


% generate the d' with the values < 95th percentile set to zero
% makes this regardless of plotting it. 
for c = 1:length(goods)
    ds = squeeze(dprime_shuffle(c,:,:))';   %for cell X, all iterationsxtime

ds_quant = quantile(ds,0.95, 1);
ds_sem = sem(ds,'omitnan');
ds_std = std(ds,'omitnan');

within = dprime(c,:)<ds_quant & dprime(c,:)>mean(ds)-ds_quant;
dprime_sig = dprime(c,:);
dprime_sig(within)= 0;

dprime_sig_cells(c,:) = dprime_sig;
end

if nargin<10
else
%keyboard
load("\\constantinoplelab.cns.nyu.edu\server\Maggie\Software\Matlab\Colormaps\TriColor.mat")
cluster = getclusternames(all_SU,goods);
figure
T = find(xvec>=win(1),1):find(xvec<=win(2),1,'last');
subplot(1,3,1)
T2 = xvec>=0 & xvec<=1;
[~,i] = sort(mean(dprime(:,T2),2),'descend');
imagesc(xvec(T),1:length(goods),dprime(i,T))
colorbar; clim([-1 1]);colormap(TriColor)
yticks(1:length(i));yticklabels(cluster(i))
title('d prime')
subplot(1,3,2)
imagesc(xvec(T),1:length(goods),dprime(:,T))
subplot(1,3,3)
imagesc(xvec(T),1:length(goods),d(i,T))
colorbar; clim([-4 4])%clim([min(d(:,T2),[],'all'),max(d(:,T2),[],'all')])
title('difference')

figure
T = find(xvec>=win(1),1):find(xvec<=win(2),1,'last');
imagesc(xvec(T),1:length(goods),dprime(:,T))
yticks(1:4:length(i));yticklabels(cluster(1:4:end))
colorbar; clim([-1 1]);colormap(TriColor)
title('not sorted dprime')


% %% individual cell
% c = 1;
% ds = squeeze(dprime_shuffle(c,:,:))';   %for cell X, all iterationsxtime
% 
% ds_quant = quantile(ds,0.95, 1);
% ds_sem = sem(ds,'omitnan');
% ds_std = std(ds,'omitnan');
% 
% % % Just shuffle comparison
% % figure
% % plot(xvec,mean(ds,'omitnan'),'r');hold on
% % shadedErrorBar(xvec,mean(ds,'omitnan'),ds_sem,'lineprops',{'color','k'})
% % shadedErrorBar(xvec,mean(ds,'omitnan'),ds_quant-mean(ds,'omitnan'),'lineprops',{'color','b'})
% % shadedErrorBar(xvec,mean(ds,'omitnan'),ds_std,'lineprops',{'color','r'})
% 
% figure
% tiledlayout(1,2)
% nexttile
% plot(xvec(T),dprime(c,T),'b')
% shadedErrorBar(xvec(T),mean(ds(:,T),'omitnan'),ds_quant(:,T)-mean(ds(:,T),'omitnan'),'lineprops',{'color','k'})
% 
% % remove values within the 95th percentile for the individual cell
% within = dprime(c,:)<ds_quant & dprime(c,:)>mean(ds)-ds_quant;
% dprime_sig = dprime(c,:);
% dprime_sig(within)= 0;
% 
% nexttile
% plot(xvec(T),dprime_sig(T))

%% all cells

figure
imagesc(xvec(T),1:length(goods),dprime_sig_cells)
colormap(TriColor);colorbar;clim([-2 2])
cluster = getclusternames(all_SU,goods);
yticks(1:length(goods));yticklabels(cluster);
box off; set(gca,'TickDir','out'); 
xlabel('Time (s)');xline(0,'--');
title({'d prime for individual clusters','Only values greater than the 95% confidence interval'})

end