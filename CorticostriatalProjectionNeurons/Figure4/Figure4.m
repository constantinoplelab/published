%% Figure 4:  OFC→CPi encoding heterogeneity supports categorical encoding of blocks from rewards

%load data
datpath = 'C:\Users\mld9131\Documents\GitHub\Papers\MLD\Figure2\Data';
load(fullfile(datpath,'OFC-CPi_projection_neurons_non-Stimulated.mat'));
CPi = 1:length(all_SU_DS);

% initialize
xvec_200 = -4:0.2:8;
xvec_100 = -4:0.1:8;
xvec = -4:0.05:8;
vcolors = [.6 .2 .9; .1 .5 .9; 0 .4 0; 1 .65 0; 1 0 0]; %need to play with colors to be color blind friendly
alignto = {'COFF', 'SON','SOFF','Rew','Opt','CON'};


%get d' of OFC-CPi projectors using 200ms bins
win = [-2,6];
trials_CPi = volbin_all_cells(all_SU_DLS_200,CPi,all_S_DS,all_index_DS,win,'Rew',0);

signedDprime = 'signed';
event = 'Rew';
shufnum = 100;
[dprime_DS,sig_DS] = vol_dprime(trials_CPi(:,[1,3]),CPi,all_S_DS,all_SU_DLS_200,all_index_DS,event,win,signedDprime,shufnum); %dprime
T = find(xvec_200>=win(1),1):find(xvec_200<=win(2),1,'last');
times = xvec_200(T);
c = isnan(dprime_DS);
d = dprime_DS;
d(c) = 0;
ind_CPi = plotdprime(d,T,xvec_200,0,1);
close()
CPi = CPi(ind_CPi);

figure('color','w')
tiledlayout(4,5)

%% A Positive d' neurons
z = 5;  %number of groups
b = 1;  %which group to study
a = (length(CPi)/z);
section  = 1+a*(b-1):a*(b-1)+a;

T = find(xvec_100>=win(1),1):find(xvec_100<=win(2),1,'last');
pvalue = CPi(section);
[avgSpk] = Vol_all_cells(all_SU_DS,pvalue,all_S_DS,all_index_DS,[-1,3],'Rew',0,1:length(pvalue),0,1);

% Firing rates
    nexttile
    for k = 4
        for r = 1:5
            plot(xvec,mean(avgSpk.(alignto{k}){r}),'color',vcolors(r,:))
            hold on
        end
    
        ylim([0 15])
        xlim([-1 3])
        set(gca,'box','off','tickdir','out')
        xlabel('Time from Reward')
        ylabel('Firing rate (Hz)')
        legend('5','10','20','40','80','Location','southeastoutside')
        title('Postive value encoding neurons')
    end

    % d' 
    [pvalued4080,psig4080] = dprimespecifyvol(pvalue,all_SU_DLS_100,all_S_DS,all_index_DS,[-1 3],1:length(pvalue),4,5,0,event);
    [pvalue1020,psig1020] = dprimespecifyvol(pvalue,all_SU_DLS_100,all_S_DS,all_index_DS,[-1 3],1:length(pvalue),1,2,0,event);

    nexttile
    shadedErrorBar(xvec_100(T),mean(pvalued4080(:,T),'omitnan'),sem(pvalued4080(:,T),'omitnan'),'LineProps',{'color',[1 0.1 0.1],'linewidth',1})
    shadedErrorBar(xvec_100(T),mean(pvalue1020(:,T),'omitnan'),sem(pvalue1020(:,T),'omitnan'),'LineProps',{'color',[0.1 0.1 1],'linewidth',1})
    hold on
    pvalued4080_2 = pvalued4080;
    pvalued4080_2(:,psig4080==0) = NaN;
    plot(xvec_100(T),mean(pvalued4080_2(:,T),'omitnan'),'color',[1 0.1 0.1],'linewidth',2)
    pvalue1020_2 = pvalue1020;
    pvalue1020_2(:,psig1020==0) = NaN;
    plot(xvec_100(T),mean(pvalue1020_2(:,T),'omitnan'),'color',[0.1 0.1 1],'linewidth',2)
    xline(0,'--')
    xlim([-0.5 1.5])
    title('Postive value encoding neurons')
    legend('40/80','5/10','Location','southeastoutside')
    set(gca,'box','off','tickdir','out')
    ylim([-.45 .4])
    yticks(-0.4:0.2:0.4)



%% B Negative d' neurons
b = 5;  %which group to study
section  = round(1+a*(b-1):a*(b-1)+a);
nvalue = CPi(section);
avgSpk_neg = Vol_all_cells(all_SU_DS,nvalue,all_S_DS,all_index_DS,[-1,3],'Rew',0,1:length(nvalue),0,1);
nexttile
for k = 4
    for r = 1:5
        plot(xvec,mean(avgSpk_neg.(alignto{k}){r}),'color',vcolors(r,:))
        hold on
    end

    ylim([0 15])
    xlim([-1 3])
    set(gca,'box','off','tickdir','out')
    xlabel('Time from Reward')
    ylabel('Firing rate (Hz)')
    legend('5','10','20','40','80','Location','southeastoutside')
    title('Negative value encoding neurons')
end

nvalued4080 = [];
nsig4080 = [];
nvalue1020 = [];
nsig1020 = [];
[nvalued4080,nsig4080] = dprimespecifyvol(nvalue,all_SU_DLS_100,all_S_DS,all_index_DS,[-1 3],1:length(nvalue),4,5,0,event);
[nvalue1020,nsig1020] = dprimespecifyvol(nvalue,all_SU_DLS_100,all_S_DS,all_index_DS,[-1 3],1:length(nvalue),1,2,0,event);

nexttile
shadedErrorBar(xvec_100(T),mean(nvalued4080(:,T),'omitnan'),sem(nvalued4080(:,T),'omitnan'),'LineProps',{'color',[1 0.1 0.1],'linewidth',1})
shadedErrorBar(xvec_100(T),mean(nvalue1020(:,T),'omitnan'),sem(nvalue1020(:,T),'omitnan'),'LineProps',{'color',[0.1 0.1 1],'linewidth',1})
hold on
nvalued4080_2 = nvalued4080;
nvalued4080_2(:,nsig4080==0) = NaN;
plot(xvec_100(T),mean(nvalued4080_2(:,T),'omitnan'),'color',[1 0.1 0.1],'linewidth',2)
nvalue1020_2 = nvalue1020;
nvalue1020_2(:,nsig1020==0) = NaN;
plot(xvec_100(T),mean(nvalue1020_2(:,T),'omitnan'),'color',[0.1 0.1 1],'linewidth',2)
xline(0,'--')
xlim([-0.5 1.5])
title('Negative value encoding neurons')
legend('40/80','5/10','Location','southeastoutside')
ylim([-.7 .3])
yticks(-0.8:0.2:0.2)

nexttile
axis off

%% C Example neuron with max firing rate trials at Offer cue and reward
nexttile
a = find(cell2mat(all_index_DS(:,3))==332); %332,318(4),338
aa = DLS(a)
S = all_S_DS{all_index_DS{a,4}};
SU = all_SU_DS{a};
V510 = S.RewardAmount<20 & S.hits==1;
V4080 = S.RewardAmount>20 & S.hits==1;
indx = all_index_DS{a,4};
[Coffmax,coffmaxi] = max(SU.hmat.COFF(:,T),[],'all');
hmat = SU.hmat.COFF(:,T);
[x,y] = find(hmat==Coffmax,1);

[rewmax,rewmaxi] = max(SU.hmat.Rew(:,T),[],'all');
hmat2 = SU.hmat.Rew(:,T);
[x2,y2] = find(hmat2==rewmax,1);

win2 = [-1 3];
T2 = find(xvec>=win2(1),1):find(xvec<=win2(2),1,'last');


plot(xvec(T2),SU.hmat.COFF(x,T2),'k')
hold on;
plot(xvec(T2),SU.hmat.Rew(x2,T2),'color',[0.5 0.5 0.5])
set(gca, 'box','off','tickdir','out');
xline(0,'--')
fill([0 0 3 3], [0 40 40 0],...
    [0.5 0.5 0.5], facealpha=0.15, edgecolor='none')
legend({'COFF','Rew'})
xlabel('Time from event')
ylabel('Firing Rate (Hz)')



win = [0 3];
T = find(xvec>=win(1),1):find(xvec<=win(2),1,'last');

mFR = [];
rFR = [];
for i = 1:length(all_SU_DS)
    mFR(i) = max(all_SU_DS{i}.hmat.COFF(:,T),[],'all');
    rFR(i) = max(all_SU_DS{i}.hmat.Rew(:,T),[],"all");
end

[h,p] = ttest(mFR,rFR); %paired ttest
[p,h] = signrank(mFR-rFR)
adtest(mFR-rFR)


nexttile 
plot(log(rFR),log(mFR),'.k')
hold on
 plot([0,5],[0,5],'color',[0.5 0.5 0.5])
ylabel('Max log(FR) during offer cue (Hz)')
xlabel('Max log(FR) during reward (Hz)')
set(gca, 'Box','off','tickdir','out')
title(strcat("signed rank, p=", string(p)))
xlim([3 5])
ylim([3 5])
axis square

nexttile
dist = histfit(log(rFR)-log(mFR),20);
dist(1).FaceColor = [0.75 0.75 0.75];
dist(1).EdgeColor = 'none';     % Remove the lines/outlines
dist(2).Color = 'k'; 
xline(mean(log(rFR) - log(mFR)), 'k')
xline(0, '--k')
set(gca,'box', 'off','TickDir','out','Color','w')
set(gcf,'color','w')
xticks([-1.5,0,1.5])
axis square
ylabel('Count')
xlabel('\Delta log(FR)')
title('inset')

%% D
% Harsha figure
nexttile
axis off
nexttile
axis off


%% E FSI waveform probability density
load('\\constantinoplelab.cns.nyu.edu\server2\PhysiologyData\EphysTable.mat');

datpath = 'C:\Users\mld9131\Documents\GitHub\Papers\MLD\Figure2\Data';
load(fullfile(datpath,'OFC-untagged_neurons_non-Stimulated.mat'),'all_index')


totalsessions = max([all_index{:,4}]);
counter = 1;
for i =1:totalsessions
    % find the sessions in the ETable
    indx = find([all_index{:,4}]==i,1);
    etable_idx = string(ETable.ratname)==all_index{indx,1} & ETable.sessiondate==all_index{indx,2} & ETable.stimulation==0 & string(ETable.recording_site)=='OFC';

    if sum(etable_idx) >1
        keyboard
    end

    % load the S struct
    load(fullfile(ETable.savepath{etable_idx},ETable.matfile{etable_idx}),'S')

    % add waveform and duration to the hmat
    if isfield(S,'templateWaveform')

        ops = make_ops(ETable(etable_idx,:));
        if ops.probe_fs==0
            keyboard
        end

        for ii = 1:length(S.templateWaveform(:,1))
            Waveform(counter,:) = S.templateWaveform(ii,:);
            Duration(counter,:) = S.templateDuration(ii,:);
            SF(counter,:) = ops.probe_fs;
            counter = counter+1;
        end    
    end
end

nexttile
DurMS = Duration./SF*1000;
[F,xi] = ksdensity(DurMS);
plot(xi, F, 'k')
set(gca,'TickDir','out','Box','off')
axis square
xlabel('Half-width duration (ms)');xlim([0 1.5]);xticks(0:0.5:1.5)
ylabel('Probability density'); yticks(0:1:4)
xline(0.4, '--k')


%% F FSI waveform comparison
nexttile
shadedErrorBar(([-41:40]./30000*1000),mean(Waveform(Duration<12,:),'omitnan'),std(Waveform(Duration<12,:),'omitnan'))
shadedErrorBar(([-41:40]./30000*1000),mean(Waveform(Duration>=12,:),'omitnan'),std(Waveform(Duration>=12,:),'omitnan'),'lineprops','r')
% shadedErrorBar(([-41:40]./30000*1000),mean(Waveform(Duration<12,:),'omitnan'),sem(Waveform(Duration<12,:),'omitnan'))
% shadedErrorBar(([-41:40]./30000*1000),mean(Waveform(Duration>=12,:),'omitnan'),sem(Waveform(Duration>=12,:),'omitnan'),'lineprops','r')

xlim([-0.5 1.5])
legend(strcat('Fast Spiker n=',string(sum(Duration<12))),strcat('Regular Spiker n= ',string(sum(Duration>12))),'Location','southeast');
set(gca,'TickDir','out','Box','off')
xlabel('Time (ms)')
ylabel('Voltage')

%% G Positive d' fast-spiking neurons
clearvars DLS trials_DS

datpath = 'C:\Users\mld9131\Documents\GitHub\Papers\MLD\Figure2\Data';
load(fullfile(datpath,'putative_FSI.mat'));
DLS = 1:length(all_SU_FSI_50);

xvec_200 = -4:0.2:8;
xvec_100 = -4:0.1:8;
xvec = -4:0.05:8;
bins = 0.2;

win = [-2,6];
trials_DS = volbin_all_cells(all_SU_FSI_200,DLS,all_S_FSI,all_index_FSI,win,'Rew',0);

signedDprime = 'signed';
event = 'Rew';
shufnum = 100;
[dprime_DS,sig_DS,DLS2] = vol_dprime(trials_DS(:,[1,3]),DLS,all_S_FSI,all_SU_FSI_200,all_index_FSI,event,win,signedDprime,shufnum); %dprime
T = find(xvec_200>=win(1),1):find(xvec_200<=win(2),1,'last');
times = xvec_200(T);
c = isnan(dprime_DS);
d = dprime_DS;
d(c) = 0;
ind_DS = plotdprime(d,T,xvec_200,0,1);
close()
DLS = DLS2(ind_DS);


z = 5;  %number of groups
b = 1;  %which group to study
a = (length(DLS)/z);
section  = 1+a*(b-1):a*(b-1)+a; 

T = find(xvec_100>=win(1),1):find(xvec_100<=win(2),1,'last');
pvalue = DLS(section);
[avgSpk] = Vol_all_cells(all_SU_FSI_50,pvalue,all_S_FSI,all_index_FSI,[-1,3],'Rew',0,1:length(pvalue),0,1);

% Firing rates
nexttile
for k = 4
    for r = 1:5
        plot(xvec,mean(avgSpk.(alignto{k}){r}),'color',vcolors(r,:))
        hold on
    end

    ylim([9 23])
    xlim([-1 3])
    set(gca,'box','off','tickdir','out')
    xlabel('Time from Reward')
    ylabel('Firing rate (Hz)')
    legend('5','10','20','40','80','Location','southeastoutside')
    title('Postive value encoding neurons')
end



% d'
[pvalued4080,psig4080] = dprimespecifyvol(pvalue,all_SU_FSI_100,all_S_FSI,all_index_FSI,[-1 3],1:length(pvalue),4,5,0,event);
[pvalue1020,psig1020] = dprimespecifyvol(pvalue,all_SU_FSI_100,all_S_FSI,all_index_FSI,[-1 3],1:length(pvalue),1,2,0,event);

nexttile
shadedErrorBar(xvec_100(T),mean(pvalued4080(:,T),'omitnan'),sem(pvalued4080(:,T),'omitnan'),'LineProps',{'color',[1 0.1 0.1],'linewidth',1})
shadedErrorBar(xvec_100(T),mean(pvalue1020(:,T),'omitnan'),sem(pvalue1020(:,T),'omitnan'),'LineProps',{'color',[0.1 0.1 1],'linewidth',1})
hold on
pvalued4080_2 = pvalued4080;
pvalued4080_2(:,psig4080==0) = NaN;
plot(xvec_100(T),mean(pvalued4080_2(:,T),'omitnan'),'color',[1 0.1 0.1],'linewidth',2) 
pvalue1020_2 = pvalue1020;
pvalue1020_2(:,psig1020==0) = NaN;
plot(xvec_100(T),mean(pvalue1020_2(:,T),'omitnan'),'color',[0.1 0.1 1],'linewidth',2)
xline(0,'--')
xlim([-0.5 1.5])
title('Positive value encoding neurons')
legend('40/80','5/10')
set(gca,'box','off','tickdir','out')
ylim([-.1 .2])
yticks(-0.4:0.2:0.4)



%% H Overlay of OFC-CPi and FSI
nexttile
load('Z:\Maggie\Papers\Physiology\Data\Rebuttle\FSI.mat')
load('Z:\Maggie\Papers\Physiology\Data\Rebuttle\DS.mat')
win = [-0.5 1.5];
T = find(xvec_100>=win(1),1):find(xvec_100<=win(2),1,'last');
b = 1;
bfsi = 1;
shadedErrorBar(xvec_100(T),mean(pvalue4080_DS{b}(:,T),'omitnan'),sem(pvalue4080_DS{b}(:,T),'omitnan'),'LineProps',{'color',[0.7 0 0],'linewidth',1})
shadedErrorBar(xvec_100(T),mean(pvalue4080_FSI{bfsi}(:,T),'omitnan'),sem(pvalue4080_FSI{bfsi}(:,T),'omitnan'),'LineProps',{'color',[1 0.1 0.1],'linewidth',1})
hold on
plot(xvec_100(T),mean(pvalue4080_sig_DS{b}(:,T),'omitnan'),'color',[0.7 0 0],'linewidth',3) 
plot(xvec_100(T),mean(pvalue4080_sig_FSI{bfsi}(:,T),'omitnan'),'color',[1 0.1 0.1],'linewidth',3) 
xline(0,'--')
xlim([-0.5, 1.5])
legend('OFC-CPi','FSI')

set(gca,'box','off','TickDir','out'); xlabel('Time (s)'); ylabel('d prime');



%% I schematic is generated in illustrator



%% J Toy data for model schematic
x = 1:5;
y = 5:5:25;
% yi = [-7.5 -2.5 0 2.5 7.5];
% yi = [0 0 2 4 6];
% yi = [0 0 2 6 10];
yi = [0 0 2 6 10];

nexttile
plot(x,y,'.k')
title('Graded inputs')
ylim([0 27])

nexttile
plot(x,yi,'.r')
title('FSI activity')

nexttile
plot(x,y-yi,'.','color',[0.5 0.5 0.5])
yline(0,'--')
title('Postive value encoding OFC-CPi activity')

for s = 16:18
    nexttile(s)
    xlim([0 6])
    xnum = {'5','10','20','40','80'};
    xticks(1:5)
    xticklabels(xnum); xlabel('Reward offer');
    ylabel('Firing rate (Hz)')
    axis square
    set(gca,'TickDir','out', 'Box','off')
end

