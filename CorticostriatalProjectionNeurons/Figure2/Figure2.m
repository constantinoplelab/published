%% Figure2
datpath = 'C:\Users\mld9131\Documents\GitHub\Papers\MLD\Figure1\Data';
load(fullfile(datpath,'ratList.mat'));

%% Optogenetics
% M: schematic of injection site and task event
load(fullfile(datpath,'CPi_opto.mat'))

% N example OFC-CPi rat
nexttile
a = 2;
rat = ConA{a}.RatName;
[hic, loc, mc] = blocks(ConA{a});
[hie, loe, me] = blocks(ExpA{a});
h = shadedErrorBar(1:5,mc.wt,mc.er,'lineprops',{'color',[0,0,0],'LineWidth',1});
set(h.edge,'LineStyle','none')
hold on
h= shadedErrorBar(1:5,me.wt,me.er,'lineprops',{'color',[0.5,0.5,0.5],'LineWidth',1});
set(h.edge,'LineStyle','none')

xlabel('Reward offer (\mu L)'); xlim([0.5 5.5]); xticks(1:5); set(gca,'TickDir','out');xticklabels({'5','10','20','40','80'});
ylim([7 12]); ylabel('Wait time (s)'); yticks(7:12)
title('Example OFC-CPi rat')


nexttile
h = shadedErrorBar(1:5,hic.wt,hic.er,'lineprops',{'color','r','LineWidth',1});
set(h.edge,'LineStyle','none')
hold on
h = shadedErrorBar(1:5,hie.wt,hie.er,'lineprops',{'color',[243/256,127/256,113/256],'LineWidth',1});
set(h.edge,'LineStyle','none')

h = shadedErrorBar(1:5,loc.wt,loc.er,'lineprops',{'color',[0,0,1],'LineWidth',1});
set(h.edge,'LineStyle','none')
h = shadedErrorBar(1:5,loe.wt,loe.er,'lineprops',{'color',[0/256,173/256,238/256],'LineWidth',1});
set(h.edge,'LineStyle','none')
xlabel('Reward offer (\mu L)'); xticks(1:5); set(gca,'TickDir','out');xticklabels({'5','10','20','40','80'});
ylim([7 12]); ylabel('Wait time (s)'); yticks(7:12)


% O population CPi rats
nexttile
[highc,highe,lowc,lowe,mixc,mixe]=wt_notsided(ExpA,ConA,Raw_data);
h = shadedErrorBar(1:5, mean(mixc,'omitnan'),sem(mixc,'omitnan'),'lineprops',{'color',[0,0,0],'LineWidth',1});
set(h.edge,'LineStyle','none')
h = shadedErrorBar(1:5, mean(mixe,'omitnan'),sem(mixc,'omitnan'),'lineprops',{'color',[0.5,0.5,0.5],'LineWidth',1});
set(h.edge,'LineStyle','none')
ylabel('Wait time (z-scored)')
xticks(1:5); xticklabels([5,10, 20, 40, 80]); xlim([0.5 5.5]); xlabel('Offered reward (uL)')
set(gca,'TickDir','out');box off;
ylim([0.84 1.32]),yticks(0.8:0.1:1.3)

for v = 1:5
[pvalue(v),z(v)] = signrank(mixe(:,v)-mixc(:,v));
end
hold on;
plot(find(z==1),z(z==1)+0.3,'*k')


nexttile
h = shadedErrorBar(1:5, mean(highc,'omitnan'),sem(highc,'omitnan'),'lineprops',{'color',[1,0,0],'LineWidth',1});
set(h.edge,'LineStyle','none')
h = shadedErrorBar(1:5, mean(highe,'omitnan'),sem(highe,'omitnan'),'lineprops',{'color',[243/256,127/256,113/256],'LineWidth',1});
set(h.edge,'LineStyle','none')
hold on;
pvalue = [];
z = [];
for v = 3:5
[pvalue(v),z(v)] = signrank(highe(:,v)-highc(:,v));
end

plot(find(z==1),z(z==1)+0.3,'*k')

h = shadedErrorBar(1:5, mean(lowc,'omitnan'),sem(lowc,'omitnan'),'lineprops',{'color',[0,0,1],'LineWidth',1});
set(h.edge,'LineStyle','none')
h = shadedErrorBar(1:5, mean(lowe,'omitnan'),sem(lowe,'omitnan'),'lineprops',{'color',[0/256,173/256,238/256],'LineWidth',1});

pvalue = [];
z = [];
for v = 1:3
[pvalue(v),z(v)] = signrank(lowe(:,v)-lowc(:,v));
end

plot(find(z==1),z(z==1)+0.3,'*k')

set(h.edge,'LineStyle','none')
ylabel('Wait time (z-scored)')
xticks(1:5); xticklabels([5,10, 20, 40, 80]); xlim([0.5 5.5]); xlabel('Offered reward (uL)')
set(gca,'TickDir','out');box off;
ylim([0.84 1.32]),yticks(0.8:0.1:1.3)


% P Delta Wait time
nexttile

lowchange = mean(lowe-lowc,2,'omitnan');
highchange = mean(highe-highc,2,'omitnan');
mixchange = mean(mixe-mixc,2,'omitnan');

plot(2.25, mean(mixchange,'omitnan'),'_k');
hold on;
errorbar(2.25,mean(mixchange,'omitnan'),sem(mixchange,'omitnan'),'-k','CapSize',0,'LineWidth',1)

plot(1.25, mean(lowchange,'omitnan'),'-b')
errorbar(1.25,mean(lowchange,'omitnan'),sem(lowchange,'omitnan'),'_b','CapSize',0,'LineWidth',1)

plot(3.25, mean(highchange,'omitnan'),'-r')
errorbar(3.25,mean(highchange,'omitnan'),sem(highchange,'omitnan'),'_r','CapSize',0,'LineWidth',1)
box off
xticks(1:3);xticklabels({'low','mix','high'}), set(gca,'TickDir','out')
yline(0,'--');xlim([0.5, 3.5]); ylabel(' \Delta Wait time (z-score)'); xlabel('Block');
ylim([-0.15 0.01]), yticks([-0.15:0.05:0.05])

% Friedman non-parametric repeated measures anova test 
% includes each volume
lowchange_v = lowe(:,1:3)-lowc(:,1:3);
highchange_v = highe(:,3:5)-highc(:,3:5);
mixchange_v = mixe(:,2:4)-mixc(:,2:4);
lowchange_v = lowchange_v';
highchange_v = highchange_v';
mixchange_v = mixchange_v';

mixchange_v = reshape(mixchange_v,1,[])';
lowchange_v = reshape(lowchange_v,1,[])';
highchange_v = reshape(highchange_v,1,[])';

T = [lowchange_v,mixchange_v,highchange_v];
%[p{4},tbl,stats] = friedman(T,3);

[p{1},h] = signrank(lowchange_v,mixchange_v);
[p{2},h] = signrank(lowchange_v,highchange_v);
[p{3},h] = signrank(mixchange_v,highchange_v);

nexttile
nexttile
nexttile

% Q  example OFC-CPr rat
% [ExpA,ConA,Raw_data] = Opto_dailycheck('VS',1,0);
% save(fullfile(datpath,'CPr_opto.mat'),'ConA','ExpA','Raw_data');
load(fullfile(datpath,'CPr_opto.mat'))

nexttile
a = 1;
rat = ConA{a}.RatName;
[hic, loc, mc] = blocks(ConA{a});
[hie, loe, me] = blocks(ExpA{a});
h = shadedErrorBar(1:5,mc.wt,mc.er,'lineprops',{'color',[0,0,0],'LineWidth',1});
set(h.edge,'LineStyle','none')
hold on
h= shadedErrorBar(1:5,me.wt,me.er,'lineprops',{'color',[0.5,0.5,0.5],'LineWidth',1});
set(h.edge,'LineStyle','none')

xlabel('Reward offer (\mu L)'); xlim([0.5 5.5]); xticks(1:5); set(gca,'TickDir','out');xticklabels({'5','10','20','40','80'});
ylim([8 16]); ylabel('Wait time (s)'); yticks(8:16)
title('Example OFC-CPr rat')


nexttile
h = shadedErrorBar(1:5,hic.wt,hic.er,'lineprops',{'color','r','LineWidth',1});
set(h.edge,'LineStyle','none')
hold on
h = shadedErrorBar(1:5,hie.wt,hie.er,'lineprops',{'color',[243/256,127/256,113/256],'LineWidth',1});
set(h.edge,'LineStyle','none')

h = shadedErrorBar(1:5,loc.wt,loc.er,'lineprops',{'color',[0,0,1],'LineWidth',1});
set(h.edge,'LineStyle','none')
h = shadedErrorBar(1:5,loe.wt,loe.er,'lineprops',{'color',[0/256,173/256,238/256],'LineWidth',1});
set(h.edge,'LineStyle','none')
xlabel('Reward offer (\mu L)');  xlim([0.5 5.5]); xticks(1:5); set(gca,'TickDir','out');xticklabels({'5','10','20','40','80'});
ylim([8 16]); ylabel('Wait time (s)'); yticks(8:16)


% R population CPr rats
nexttile
[highc,highe,lowc,lowe,mixc,mixe]=wt_notsided(ExpA,ConA,Raw_data);
h = shadedErrorBar(1:5, mean(mixc,'omitnan'),sem(mixc,'omitnan'),'lineprops',{'color',[0,0,0],'LineWidth',1});
set(h.edge,'LineStyle','none')
h = shadedErrorBar(1:5, mean(mixe,'omitnan'),sem(mixc,'omitnan'),'lineprops',{'color',[0.5,0.5,0.5],'LineWidth',1});
set(h.edge,'LineStyle','none')
ylabel('Wait time (z-scored)')
xticks(1:5); xticklabels([5,10, 20, 40, 80]); xlim([0.5 5.5]); xlabel('Offered reward (uL)')
set(gca,'TickDir','out');box off;
ylim([0.84 1.32]),yticks(0.8:0.1:1.3)
title("OFC-CPr rats n = 7")

for v = 1:5
[pvalue(v),z(v)] = signrank(mixe(:,v)-mixc(:,v));
end
hold on;
plot(find(z==1),z(z==1)+0.3,'*k')


nexttile
h = shadedErrorBar(1:5, mean(highc,'omitnan'),sem(highc,'omitnan'),'lineprops',{'color',[1,0,0],'LineWidth',1});
set(h.edge,'LineStyle','none')
h = shadedErrorBar(1:5, mean(highe,'omitnan'),sem(highe,'omitnan'),'lineprops',{'color',[243/256,127/256,113/256],'LineWidth',1});
set(h.edge,'LineStyle','none')
hold on;
pvalue = [];
z = [];
for v = 3:5
[pvalue(v),z(v)] = signrank(highe(:,v)-highc(:,v));
end

plot(find(z==1),z(z==1)+0.3,'*k')

h = shadedErrorBar(1:5, mean(lowc,'omitnan'),sem(lowc,'omitnan'),'lineprops',{'color',[0,0,1],'LineWidth',1});
set(h.edge,'LineStyle','none')
h = shadedErrorBar(1:5, mean(lowe,'omitnan'),sem(lowe,'omitnan'),'lineprops',{'color',[0/256,173/256,238/256],'LineWidth',1});

pvalue = [];
z = [];
for v = 1:3
[pvalue(v),z(v)] = signrank(lowe(:,v)-lowc(:,v));
end

plot(find(z==1),z(z==1)+0.3,'*k')

set(h.edge,'LineStyle','none')
ylabel('Wait time (z-scored)')
xticks(1:5); xticklabels([5,10, 20, 40, 80]); xlim([0.5 5.5]); xlabel('Reward offer (\mu L)');
set(gca,'TickDir','out');box off;
ylim([0.84 1.32]),yticks(0.8:0.1:1.3)


% S Delta Wait time
nexttile

lowchange = mean(lowe-lowc,2,'omitnan');
highchange = mean(highe-highc,2,'omitnan');
mixchange = mean(mixe-mixc,2,'omitnan');

plot(2.25, mean(mixchange,'omitnan'),'_k');
hold on;
errorbar(2.25,mean(mixchange,'omitnan'),sem(mixchange,'omitnan'),'-k','CapSize',0,'LineWidth',1)

plot(1.25, mean(lowchange,'omitnan'),'-b')
errorbar(1.25,mean(lowchange,'omitnan'),sem(lowchange,'omitnan'),'_b','CapSize',0,'LineWidth',1)

plot(3.25, mean(highchange,'omitnan'),'-r')
errorbar(3.25,mean(highchange,'omitnan'),sem(highchange,'omitnan'),'_r','CapSize',0,'LineWidth',1)
box off
xticks(1:3);xticklabels({'Low','Mix','High'}), set(gca,'TickDir','out')
yline(0,'--');xlim([0.5, 3.5]); ylabel(' \Delta Wait time (z-score)'); xlabel('Block');
ylim([-0.15 0.01]), yticks([-0.15:0.05:0.05])

% Friedman non-parametric repeated measures anova test 
% includes each volume
lowchange_v = lowe(:,1:3)-lowc(:,1:3);
highchange_v = highe(:,3:5)-highc(:,3:5);
mixchange_v = mixe(:,2:4)-mixc(:,2:4);
lowchange_v = lowchange_v';
highchange_v = highchange_v';
mixchange_v = mixchange_v';

mixchange_v = reshape(mixchange_v,1,[])';
lowchange_v = reshape(lowchange_v,1,[])';
highchange_v = reshape(highchange_v,1,[])';

T = [lowchange_v,mixchange_v,highchange_v];
%[p{4},tbl,stats] = friedman(T,3);

[p{1},h] = signrank(lowchange_v,mixchange_v);
[p{2},h] = signrank(lowchange_v,highchange_v);
[p{3},h] = signrank(mixchange_v,highchange_v);

% Not shown but included in text
% Compare Delta WT between CPi and CPr stimulated animals
p = compareDeltaWT_CPiCPr;


nexttile
nexttile
nexttile

%% Biased Model
% T hand drawn schematic of biased priors

Nsims = 100;
load(fullfile(datpath,'Bias_Model_data_new.mat'))

clearvars WT

[WT{1}, WT{2},WT{3}] = blocks(A_mdl);
[WT_bias{1}, WT_bias{2},WT_bias{3}] = blocks(A_bias);
c = {'r', 'b','k'};
c_bias = [[243/256,127/256,113/256]; [0/256,173/256,238/256]; 0.5,0.5,0.5];

% U model WT
nexttile
for i = 1:3
    h = shadedErrorBar(1:5,WT{i}.wt,WT{i}.er,'lineProps',{"Color",c{i}}); hold on;
    set(h.edge,'LineStyle','none')
    h = shadedErrorBar(1:5,WT_bias{i}.wt,WT_bias{i}.er,'lineProps',{"Color",c_bias(i,:)}); hold on;
    set(h.edge,'LineStyle','none')
end
xticks(1:5); xticklabels([5,10, 20, 40, 80]); xlim([0.5 5.5]); xlabel('Offered reward (uL)')
set(gca, 'TickDir', 'out');
ylabel('Wait time (AU)'); yticklabels([]);
title(strcat('belief bias:', {' '},string(bias)))
legend({'control','opto', 'control', 'opto'},'location','southeast');

% V model transition dynamics
nexttile
% mix to low
% subplot(2,8,4)
h = shadedErrorBar(xvec,mean(mtol),std(mtol),'lineprops',{'color','b'}); hold on;set(h.edge,'LineStyle','none')
h = shadedErrorBar(xvec,mean(mtolb),std(mtolb),'lineprops',{'color',c_bias(2,:)});set(h.edge,'LineStyle','none')
title({'Wait-time', 'Mixed -> Adapt'})
set(gca,'TickDir','out');xline(0,'--')
xlim([-15 30]),xticks(-30:15:30)

% mix to high
% subplot(2,8,12)
h = shadedErrorBar(xvec,mean(mtoh),std(mtoh),'lineprops',{'color','r'});set(h.edge,'LineStyle','none')
h = shadedErrorBar(xvec,mean(mtohb),std(mtohb),'lineprops',{'color',c_bias(1,:)});set(h.edge,'LineStyle','none')
xlabel('Trials from block switch'), xline(0,'--');
set(gca,'TickDir','out');
xlim([-15 30]),xticks(-30:15:30)

%% W Animal block transition
load(fullfile(datpath,'CPi_opto.mat'))
clearvars WT WTe deltaWT deltaWTe

twin = 40;
smooth = 5;
usecausal = 1;
binSize = 2;
for a = 1:length(ConA)
    [WT.ltom(a,:), WT.htom(a,:), WT.mtol(a,:), WT.mtoh(a,:),WT.ltom_incong(a,:), WT.htom_incong(a,:)] =block_dynamics_wt_binTrials(ConA{a}, twin, binSize, smooth);
    [WTe.ltom(a,:), WTe.htom(a,:), WTe.mtol(a,:), WTe.mtoh(a,:),WTe.ltom_inconge(a,:),WTe.htom_inconge(a,:)] =block_dynamics_wt_binTrials(ExpA{a}, twin, binSize, smooth);
end
xvec = -twin:binSize:twin-binSize;

Before_lower = -10; Before_upper = 0;
Early_lower = 0; Early_upper = 5;
Early2_lower = 5; Early2_upper = 10;
Late_lower = 10; Late_upper = 20;
Later_lower = 20; Later_upper = 30;


nexttile
yl = [-0.1 0.25];
fill([0 1.6 1.6 0], [yl(1) yl(1) yl(2) yl(2)],...
    'k', facealpha=0.15, edgecolor='none')
hold on

win = xvec>=Before_lower & xvec<Before_upper;
deltaWT(:,1) = mean(WT.mtol(:,win),2,'omitnan');
deltaWTe(:,1) = mean(WTe.mtol(:,win),2,'omitnan');
errorbar(1,mean(mean(WT.mtol(:,win),2,'omitnan')),sem(mean(WT.mtol(:,win),2,'omitnan')),'_b','CapSize',0,'LineWidth',1)
hold on
errorbar(1.25,mean(mean(WTe.mtol(:,win),2,'omitnan')),sem(mean(WTe.mtol(:,win),2,'omitnan')),'_','color', c_bias(2,:),'CapSize',0,'LineWidth',1)

win = xvec>=Early_lower & xvec<Early_upper;
deltaWT(:,2) = mean(WT.mtol(:,win),2,'omitnan');
deltaWTe(:,2) = mean(WTe.mtol(:,win),2,'omitnan');
errorbar(2,mean(mean(WT.mtol(:,win),2,'omitnan')),sem(mean(WT.mtol(:,win),2,'omitnan')),'_b','CapSize',0,'LineWidth',1)
hold on
errorbar(2.25,mean(mean(WTe.mtol(:,win),2,'omitnan')),sem(mean(WTe.mtol(:,win),2,'omitnan')),'_','color', c_bias(2,:),'CapSize',0,'LineWidth',1)

win = xvec>=Early2_lower & xvec<Early2_upper;
deltaWT(:,3) = mean(WT.mtol(:,win),2,'omitnan');
deltaWTe(:,3) = mean(WTe.mtol(:,win),2,'omitnan');
errorbar(3,mean(mean(WT.mtol(:,win),2,'omitnan')),sem(mean(WT.mtol(:,win),2,'omitnan')),'_b','CapSize',0,'LineWidth',1)
hold on
errorbar(3.25,mean(mean(WTe.mtol(:,win),2,'omitnan')),sem(mean(WTe.mtol(:,win),2,'omitnan')),'_','color', c_bias(2,:),'CapSize',0,'LineWidth',1)

win = xvec>=Late_lower & xvec<Late_upper;
deltaWT(:,4) = mean(WT.mtol(:,win),2,'omitnan');
deltaWTe(:,4) = mean(WTe.mtol(:,win),2,'omitnan');
errorbar(4,mean(mean(WT.mtol(:,win),2,'omitnan')),sem(mean(WT.mtol(:,win),2,'omitnan')),'_b','CapSize',0,'LineWidth',1)
hold on
errorbar(4.25,mean(mean(WTe.mtol(:,win),2,'omitnan')),sem(mean(WTe.mtol(:,win),2,'omitnan')),'_','color', c_bias(2,:),'CapSize',0,'LineWidth',1)

win = xvec>=Later_lower & xvec<Later_upper;
deltaWT(:,5) = mean(WT.mtol(:,win),2,'omitnan');
deltaWTe(:,5) = mean(WTe.mtol(:,win),2,'omitnan');
errorbar(5,mean(mean(WT.mtol(:,win),2,'omitnan')),sem(mean(WT.mtol(:,win),2,'omitnan')),'_b','CapSize',0,'LineWidth',1)
hold on
errorbar(5.25,mean(mean(WTe.mtol(:,win),2,'omitnan')),sem(mean(WTe.mtol(:,win),2,'omitnan')),'_','color',c_bias(2,:),'CapSize',0,'LineWidth',1)


xlim([0.5 5.5]);xlabel('Trials from block switch'),xticks(1:5),xticklabels({'-10-0','0-5','5-10','10-20','20-30'})
ylabel('Wait time (z-score)')
box off, set(gca,'TickDir','out')
title('OFC-CPi rats')


for i = 1:5
[p(i),z(i)] = signrank(deltaWT(:,i),deltaWTe(:,i),'tail','right');
end
xvec = 1:5;
plot(xvec(z>0),z(z>0)*0.25,'*k')


% % Friedman non-parametric repeated measures anova test 
% % includes each volume
% T = [deltaWT;deltaWTe];
% dWT = reshape(deltaWT,55,1);
% dWTe = reshape(deltaWTe,55,1);
% T = [dWT,dWTe];
% 
% [p,tbl,stats] = friedman(T,11);
% multcompare(stats)

legend('control','opto','Location','southeast')