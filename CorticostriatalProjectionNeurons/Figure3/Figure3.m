%% Figure 3 Neuropixels
% A schematic
% B waveforms
% C rasters



% load your cells and full OFC.
    % these take a while to load and might crash computers
datpath = 'C:\Users\mld9131\Documents\GitHub\Papers\MLD\Figure2\Data';
load(fullfile(datpath,'OFC-CPi_projection_neurons_non-Stimulated.mat'));
load(fullfile(datpath,'OFC-CPr_projection_neurons_non-Stimulated.mat'));
load(fullfile(datpath,'OFC-untagged_neurons_non-Stimulated.mat'))

load(fullfile(datpath,'TriColor.mat'))

%initialize window
win = [-1,3];
xvec = -4:0.05:8;
alignto = {'COFF', 'SON','SOFF','Rew','Opt','CON'};

% initialize colors
vcolors = [.6 .2 .9; .1 .5 .9; 0 .4 0; 1 .65 0; 1 0 0]; 
Label_vol = {'5 \muL','10 \muL','20 \muL','40 \muL','80 \muL'};
CPr_color = [11/256, 187/256, 189/256];
CPi_color = [164/256, 38/256, 113/256];
load(fullfile(datpath,'TriColor.mat'))

%% Create the figure
figure('color','w')

% D reward volume PSTH
avgSpk = Vol_all_cells(all_SU_DS,1:length(all_SU_DS),all_S_DS,all_index_DS,[-1,3],'Rew',0,1,1,1);
T = find(xvec>=win(1),1):find(xvec<=win(2),1,'last');
limit = 0;

for k = 1:4
    subplot(3,5,k)
    for vol = 1:5
        plot(xvec(T),mean(avgSpk.(alignto{k}){vol}(:,T),'omitnan'),'Color',vcolors(vol,:))
        limit = [limit,mean(avgSpk.(alignto{k}){vol}(:,T),'omitnan')];
        hold on;
    end
end
subplot(3,5,1);
    ylabel('Firing rate (Hz)')
for k = 1:4
    subplot(3,5,k)
    top = max(limit(1,2:end));
    bottom = min(limit(1,2:end));
    ylim([bottom-abs(bottom/10) top+top/50])
    title(alignto{k})
    box off; xlabel('Time from event (s)');set(gca,'TickDir','out');xline(0,'--');
end
legend(Label_vol,'location','northeast')


% E reward volume PSTH CPr
avgSpk = Vol_all_cells(all_SU_VS,1:length(all_SU_VS),all_S_VS,all_index_VS,[-1,3],'Rew',0,1,1,1);
T = find(xvec>=win(1),1):find(xvec<=win(2),1,'last');
limit = 0;

for k = 1:4
    subplot(3,5,k+5)
    for vol = 1:5
        plot(xvec(T),mean(avgSpk.(alignto{k}){vol}(:,T),'omitnan'),'Color',vcolors(vol,:))
        limit = [limit,mean(avgSpk.(alignto{k}){vol}(:,T),'omitnan')];
        hold on;
    end
end
subplot(3,5,1+5);
    ylabel('Firing rate (Hz)')
for k = 1:4
    subplot(3,5,k+5)
    top = max(limit(1,2:end));
    bottom = min(limit(1,2:end));
    ylim([bottom-abs(bottom/10) top+top/50])
    title(alignto{k})
    box off; xlabel('Time from event (s)');set(gca,'TickDir','out');xline(0,'--');
    xlim([-1 3])
end
legend(Label_vol,'location','northeast')


% Volume d'
DLS = 1:length(all_SU_DS);
VS = 1:length(all_SU_VS);

win = [-2,6];
signedDprime = 'unsigned';
shufnum = 100;

%calculate volume trials
trials_DS = volbin_all_cells(all_SU_DLS_200,DLS,all_S_DS,all_index_DS,win,[],0);
trials_VS = volbin_all_cells(all_SU_VS_200,VS,all_S_VS,all_index_VS,win,[],0);
for k = 1:5
event = alignto{k};
[dprime_DS.(event),sig_DS.(event)] = vol_dprime(trials_DS(:,[1,3]),DLS,all_S_DS,all_SU_DLS_200,all_index_DS,event,win,signedDprime,shufnum); %dprime
[dprime_VS.(event),sig_VS.(event)] = vol_dprime(trials_VS(:,[1,3]),VS,all_S_VS,all_SU_VS_200,all_index_VS,event,win,signedDprime,shufnum); %dprime
T = find(xvec_200>=win(1),1):find(xvec_200<=win(2),1,'last');
times = xvec_200(T);
end

%for iterations of randomly selected cells.
trials_random = volbin_all_cells(all_SU_random_200,cells,all_S,all_index,win,[],0);
for k = 1:4
    event = alignto{k};
    [dprime_random.(event),sig_random.(event)] = vol_dprime(trials_random(:,[1,3]),cells,all_S,all_SU_random_200,all_index,event,win,signedDprime,shufnum); %dprime
    dprime_SEM.(event)=sem(dprime_random.(event),'omitnan');
end
dprime_mean = dprime_random;

for k = 1:4
    event = alignto{k};
% permutation test for significant difference between DS and the random OFC neurons
for i = 1:length(T)
pval.(event).p(i) = permute_test(dprime_mean.(event)(:,T(i)),dprime_DS.(event)(:,T(i)),100);
pval.(event).p2(i) = permute_test(dprime_mean.(event)(:,T(i)),dprime_VS.(event)(:,T(i)),100);
end
end

% Plot OFC-CPi,OFC-CPr and OFC neurons d' 5/10 v 40/80
for k = 1:4
event = alignto{k};
subplot(3,5,k+10)
shadedErrorBar(times,mean(dprime_DS.(event)(:,T),'omitnan'),sem(dprime_DS.(event)(:,T),'omitnan'),'LineProps',{'color',CPi_color,'linewidth',1.5})
shadedErrorBar(times,mean(dprime_VS.(event)(:,T),'omitnan'),sem(dprime_VS.(event)(:,T),'omitnan'),'LineProps',{'color',CPr_color,'linewidth',1.5})
shadedErrorBar(times,mean(dprime_mean.(event)(:,T),1,'omitnan'),dprime_SEM.(event)(:,T),'LineProps',{'color','k','linewidth',1.5})

% bold significant times
if sum(sig_DS.(event)(T))>0
    ds_sig = mean(dprime_DS.(event),'omitnan');
    ds_sig(sig_DS.(event)==0)=NaN;
    hold on;
    plot(times, ds_sig(T),'linewidth',3,'color',CPi_color)
    sig= pval.(event).p<0.05 & (~ismissing(ds_sig(T)));
    h = ones(1,length(T));
    plot(times(sig)+0.05,h(sig)-0.25,'.','color',CPi_color)
else
        sig = zeros(1,length(T));
end
if sum(sig_VS.(event)(T))>0
    vs_sig = mean(dprime_VS.(event),'omitnan');
    vs_sig(sig_VS.(event)==0)=NaN;
    hold on;
    plot(times, vs_sig(T),'linewidth',3,'color',CPr_color)
    sig2 = pval.(event).p2<0.05 & (~ismissing(vs_sig(T)));
    hold on
    h = ones(1,length(T));
    plot(times(sig2)+0.05,h(sig2)-0.3,'.','color',CPr_color)

else
    sig2 = zeros(1,length(T));
end
if sum(sig_random.(event)(T))>0
    random_sig = mean(dprime_mean.(event),'omitnan');
    random_sig(sig_random.(event)==0)=NaN;
    hold on
    plot(times, random_sig(T),'linewidth',3,'color','k');
end

set(gca,'Tickdir','out');


ylim([-0.25 0.8]);
xline(0,'--')
xlim([-1 3])
box off
end
subplot(3,5,11)
ylabel("Volume d'")


% Signed d' hmat
%Plot hmat of signed d' for all CPi cells
subplot(3,5,5)
[dprime_signed_DS.(event),sig_DS.(event)] = vol_dprime(trials_DS(:,[1,3]),DLS,all_S,all_SU_DLS_200,all_index,event,win,'signed',shufnum); %dprime
c = isnan(dprime_signed_DS.(event));
d = dprime_signed_DS.(event);
d(c) = 0;

v0 = d(:,xvec_200>=0 & xvec_200<=1);
v1 = find(v0==-Inf | v0==Inf);
v0(v1) = NaN;
v0 = mean(v0',1,'omitnan')';
[~,ind]=sort(v0,'descend');
imagesc(xvec_200(T),1:length(d(:,1)),d(ind,T))
xline(0, '--');
box off; set(gca,'TickDir','out');
title Reward; 
xlabel('Time from event (s)'), xlim([-1 3])
ylabel('OFC-CPi cells'); yticks([]);
clim([-1.5 1.5]); colormap(TriColor), colorbar('eastoutside','ticks',-1.5:.75:1.5,'Box','off')
set(gcf,'position',[1950,50,1850,950])


% find greater d' than 0
v0 = d(:,xvec_200>=0 & xvec_200<=1);
v1 = find(v0==-Inf | v0==Inf);
v0(v1) = NaN;
% v0 = mean(v0',1,'omitnan')'; % 53% >0 averaged across whole time
% v0 = sort(v0);
for i = 1:6
percentPositive(i) = sum(v0(:,i)>0)/sum(v0(:,i)~=0)
end
ttest(percentPositive,0.5)
mean(percentPositive)
% look for cells with significant d'


% Signed d' hmat
%Plot hmat of signed d' for all cells
subplot(3,5,10)
[dprime_signed_VS.(event),sig_VS.(event)] = vol_dprime(trials_VS(:,[1,3]),VS,all_S,all_SU_VS_200,all_index,event,win,'signed',shufnum); %dprime
c = isnan(dprime_signed_VS.(event));
d = dprime_signed_VS.(event);
d(c) = 0;

v0 = d(:,xvec_200>=0 & xvec_200<=1);
v1 = find(v0==-Inf | v0==Inf);
v0(v1) = NaN;
v0 = mean(v0',1,'omitnan')';
[~,ind]=sort(v0,'descend');
imagesc(xvec_200(T),1:length(d(:,1)),d(ind,T))
xline(0, '--');
box off; set(gca,'TickDir','out');
title Reward; 
xlabel('Time from event (s)'), xlim([-1 3])
ylabel('OFC-CPi cells'); yticks([]);
clim([-1.5 1.5]); colormap(TriColor), colorbar('eastoutside','ticks',-1.5:.75:1.5,'Box','off')
set(gcf,'position',[1950,50,1850,950])
