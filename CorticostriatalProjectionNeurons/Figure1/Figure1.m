%% Load all rat data
datpath = 'C:\Users\mld9131\Documents\GitHub\Papers\MLD\Figure1\Data';
load(fullfile(datpath,'ratList.mat'));
%datapath = 'path to the S structs'; %path to the Sstructs;

WT = struct;
[WT.high, WT.low, WT.mixed]  = deal(nan(length(ratList), 5));

% Process each rat's data
for rr = 1:length(ratList)
    fprintf('%d out of %d\n', rr, length(ratList));

    % Load Data
    fname = strcat(['ratTrial_', ratList{rr}, '.mat']);
    a = load([datapath filesep fname]);
    A = a.A;

    %detrend wait times
    A = detrendwt(A);
    
    Az = A;
    Az.wait_time(logical(A.catch)) = (A.wait_time(logical(A.catch))-...
    mean(A.wait_time(logical(A.catch) & A.block==1 & A.reward<32 & A.reward>10), 'omitnan'))./...
    std(A.wait_time(logical(A.catch)& A.block==1 & A.reward<32 & A.reward>10), 'omitnan');

    % Average wait time and SEM as a function of reward in each block
    try
    [high, low, mix, ps(rr,1)] = wtcurves(A);
    [highz, lowz, mixz, psz(rr,1)] = wtcurves(Az);
    catch
        high.wt=NaN(1,5); low.wt=NaN(1,5); mix.wt=NaN(1,5); ps(rr,1)=NaN;
        fprintf([ratList{rr}, ' does not have blocks'])
        highz.wt=NaN(1,5); lowz.wt=NaN(1,5); mixz.wt=NaN(1,5); psz(rr,1)=NaN;
    end

    WT.high(rr,:) = high.wt;
    WT.low(rr,:) = low.wt;
    WT.mixed(rr,:) = mix.wt;

    WTz.high(rr,:) = highz.wt;
    WTz.low(rr,:) = lowz.wt;
    WTz.mixed(rr,:) = mixz.wt;
    
    
 % Number of sessions and trials
    ntrials(rr,:) = [length(A.date) length(A.reward)];

% Block dynamics 
twin = 40; %trial window for wait time dynamics plot
smoothfactor = 5;
binSize = 1;

[ltom(rr,:), htom(rr,:), mtol(rr,:), mtoh(rr,:), ~, ~, ~] =...
    block_dynamics_wt_binTrials(A, twin, binSize, smoothfactor);

[expert.postLow(rr,:), expert.postHigh(rr,:), expert.postLow_q1(rr,:),  ...
    expert.postHigh_q1(rr, :), expert.postLow_q1Con(rr,:), ...
    expert.postHigh_q1Con(rr,:)] = quartileAnalysis_SS(A);

[expertz.postLow(rr,:), expertz.postHigh(rr,:), expertz.postLow_q1(rr,:),  ...
    expertz.postHigh_q1(rr, :), expertz.postLow_q1Con(rr,:), ...
    expertz.postHigh_q1Con(rr,:)] = quartileAnalysis_SS(Az);

end

%% Start the Figure
figure("Color",'white', 'PaperUnits','centimeters','PaperSize',[0.18 0.17])
tiledlayout(5,8)
xs = 1:5;

%% Figure 1 Inferential behavior
% A & B are schematics
% C. Population WT curves
pwt = signrank(WT.high(:,3), WT.low(:,3));
efs = meanEffectSize(WT.high(:,3), WT.low(:,3),Paired=true);
wt_ratio = WT.high(:,3)./WT.low(:,3);

nexttile
l(4) = shadedErrorBar(xs, mean(WT.mixed,'omitnan'),...
    sem(WT.mixed,'omitnan'));hold on;set(l(4).edge,'LineStyle','none')

l(5) = shadedErrorBar(xs, mean(WT.high,'omitnan'),...
    sem(WT.high,'omitnan'),...
    LineProps={'r'});hold on;set(l(5).edge,'LineStyle','none')

l(6) = shadedErrorBar(xs, mean(WT.low,'omitnan'),...
    sem(WT.low,'omitnan'),...
    LineProps={'b'});hold on;set(l(5).edge,'LineStyle','none')

xlim([0.5 5.5])
ylim([9 13.6]); yticks(10:1:13)

xticks(1:5)
xticklabels({'5', '10', '20', '40', '80'})

xlabel('Reward offer (μL)')
ylabel('Wait time (s)')
title(['n = ' num2str(length(ratList)),' rats'])
subtitle(strcat('p=',string(pwt)))
set(gca, 'TickDir', 'out'); box off;
for ii = 4:6
    set(l(ii).edge, 'linestyle','none')
end

% D Wait time dynamics
all_ltom = mean(ltom,'omitnan');
all_htom= mean(htom,'omitnan');
all_mtol= mean(mtol,'omitnan');
all_mtoh= mean(mtoh,'omitnan');

all_ltom_sem = sem(ltom,'omitnan');
all_htom_sem= sem(htom,'omitnan');
all_mtol_sem= sem(mtol,'omitnan');
all_mtoh_sem= sem(mtoh,'omitnan');

xvec = -twin:twin-1;

nexttile
fill([-30 -30 0 0], [-0.25 0.25 0.25 -0.25],...
    'k', facealpha=0.15, edgecolor='none')
h = shadedErrorBar(xvec, all_mtol, all_mtol_sem, 'lineprops', ...
    {'-b', 'linewidth', 1});hold on;set(h.edge,'LineStyle','none')
hold on
h = shadedErrorBar(xvec, all_mtoh, all_mtoh_sem, ...
    'lineprops', {'-r', 'linewidth', 1});hold on;set(h.edge,'LineStyle','none')
set(gca, 'TickDir', 'out'); box off;

yticks([-0.2, 0, 0.2])
ylim([-0.25 0.25])
xlabel('Trials from block switch');
ylabel('Wait time (z-scored)');
title('Mixed to Adapt')
xlim([-25 40]);
xticks(-20:20:20)


nexttile
fill([0 0 25 25], [-0.25 0.25 0.25 -0.25],...
    'k', facealpha=0.15, edgecolor='none')
h = shadedErrorBar(xvec, all_ltom, all_ltom_sem, 'lineprops', ...
    {'-b', 'linewidth', 1});hold on;set(h.edge,'LineStyle','none')
hold on
h = shadedErrorBar(xvec, all_htom, all_htom_sem, ...
    'lineprops', {'-r', 'linewidth', 1});hold on;set(h.edge,'LineStyle','none')
set(gca, 'TickDir', 'out'); box off;
yticks([-0.2, 0, 0.2])
ylim([-0.25 0.25])
xlabel('Trial from block switch');
title('Adapt to Mixed')
xlim([-25 40]);
xticks(-20:20:20)

% E Oppertunity cost schematic
nexttile
x = 0:10;
y = (2*exp(-0.2*x) + 0.1*x)-1.2;
y(7:end) = 0;
plot(x,y,'color',[1 0.45 0.25],'LineWidth',1.25)

hold on
yline(0.35,'--b')
yline(0.175,'--r')
y4 = 0:0.05:0.35;
x4 = ones(length(y4),1)*2;
plot(x4, y4, '--b')

y5 = 0:0.05:0.175;
x5 = ones(length(y5),1)*3.25;
plot(x5, y5, '--r')

ylim([-0.01 0.8])
box off
yticks([])
xticks([])
xlabel('Time')
ylabel('Value')
% E probability
nexttile
yl = [0.0 1.0];

fill([0 40 40 0], [yl(1) yl(1) yl(2) yl(2)],...
    'k', facealpha=0.15, edgecolor='none')
hold on
fill([40 80 80 40], [yl(1) yl(1) yl(2) yl(2)],...
    'b', facealpha=0.15, edgecolor='none')
fill([80 120 120 80], [yl(1) yl(1) yl(2) yl(2)],...
    'k', facealpha=0.15, edgecolor='none')
fill([120 160 160 120], [yl(1) yl(1) yl(2) yl(2)],...
    'r', facealpha=0.15, edgecolor='none')
fill([160 200 200 160], [yl(1) yl(1) yl(2) yl(2)],...
    'k', facealpha=0.15, edgecolor='none')

rng(30)
alpha0 = 0.5;
D = 20;

r = [randsample(1:5, 40, true),...
    randsample(1:3, 40, true),...
    randsample(1:5, 40, true),...
    randsample(3:5, 40, true),...
    randsample(1:5, 40, true)]';

[~, ~, V_DB1, ~, ~, G_DB1, belief] =...
    generate_deltabelief_mdl(r, alpha0, D);

plot(belief(1,:), 'k')
hold on
plot(belief(2,:), 'r')
plot(belief(3,:), 'b')
hold off
yticks([])
ylabel('P(B_t|R_t)')
title('P(B_t|R_t) ~ P(R_t|B_t)P(B_t)')

%% Inferential model
Nsims = 100;
rng(100) % Set random seed for reproducibility 

kappa_mi = 1; %.5;  
kappa_hi = 1.5; %1.2; %0.7; 
kappa_lo = 0.5;%0.8; %.3; 

D = 0.5;          %2
lambda =1;     %1
noise = 0; %.1;%2
params = [0.25 0.30 0.2 .13]; %shannon data

% Preallocate matricies
[WTbin.ltom, WTbin.htom, WTbin.mtol, WTbin.mtoh,] = deal(nan(Nsims, length(xvec)));
for  rr = 1:Nsims

    fname = strcat(['ratTrial_', ratList{rr}, '.mat']);
    a = load([datapath filesep fname]);
    A = a.A;

    A.trainingstage(1:length(A.block),1) = 9;
    A_mdl = A;
   fprintf('%d out of %d\n', rr, Nsims);

% run model
[~, A_mdl.wait_time, BlkInf, Belief, Prior, Kappa] = GenerateSynthData_Bayes(params, A,'logn', true,noise);


%repeat with noise
[WTbin.ltom(rr,:), WTbin.htom(rr,:), WTbin.mtol(rr,:), WTbin.mtoh(rr,:)] =block_dynamics_wt_binTrials(A_mdl, twin, binSize, smoothfactor);

    [infModel.postLow(rr,:), infModel.postHigh(rr,:), infModel.postLow_q1(rr,:),  ...
        infModel.postHigh_q1(rr, :), infModel.postLow_q1Con(rr,:), ...
        infModel.postHigh_q1Con(rr,:)] = quartileAnalysis_SS(A_mdl);
end

[WT_mld{1}, WT_mld{2},WT_mld{3}] = blocks(A_mdl);
c = {'r', 'b','k'};
x2 = [5 10 20 40 80];
x = log(log(x2));

% F Model wait times
nexttile
for i = 1:3
    h = shadedErrorBar(xs,WT_mld{i}.wt,WT_mld{i}.er,'lineProps',{"Color",c{i}}); hold on;
    set(h.edge,'LineStyle','none')
end

xlim([0.5 5.5])
yticks(8:13)
ylim([8.8 13])
xticks(xs);
xticklabels({'5', '10', '20', '40', '80'})
xlabel('Reward offer (μL)')
ylabel('Wait time (s)')
set(gca, 'TickDir', 'out'); box off;

% G Model Wait time dynamics
nexttile
fill([-25 -25 0 0], [-1.5 1.5 1.5 -1.5],...
    'k', facealpha=0.15, edgecolor='none')
h = shadedErrorBar(xvec,mean(WTbin.mtol),std(WTbin.mtol),'lineprops',{'color','b'}); hold on;set(h.edge,'LineStyle','none')
h= shadedErrorBar(xvec,mean(WTbin.mtoh),std(WTbin.mtoh),'lineprops',{'color','r'}); hold on;set(h.edge,'LineStyle','none')
title('Mixed to Adapt')
xlabel('Trials from block switch');
ylabel('Wait time (z-scored)');
set(gca,'TickDir','out');
yticks([-1, -0.5, 0, 0.5, 1])
ylim([-1.5 1.5])
xlim([-25 40]);
xticks(-20:20:20)
box off;

nexttile
fill([0 0 40 40], [-1.5 1.5 1.5 -1.5],...
    'k', facealpha=0.15, edgecolor='none')
h = shadedErrorBar(xvec,mean(WTbin.ltom),std(WTbin.ltom),'lineprops',{'color','b'}); hold on;set(h.edge,'LineStyle','none')
h= shadedErrorBar(xvec,mean(WTbin.htom),std(WTbin.htom),'lineprops',{'color','r'}); hold on;set(h.edge,'LineStyle','none')
title('Adapt to Mixed')
xlabel('Trials from block switch');
set(gca,'TickDir','out');
yticks([-1, -0.5, 0, 0.5, 1])
ylim([-1.5 1.5])
xlim([-25 40]);
xticks(-20:20:20)
box off;

%% Anatomy
% H-J Generated manually in illustrator

% K CPr anatomy
load(fullfile(datpath,'DCS.mat'))
CPi = CPi_data;

nexttile
imagesc(CPi')
c = colorbar;
clim([0 20])
title('OFC - CPi')
set(gca,'XColor', 'none','YColor','none','box','off','Colormap', cmap3,'TickDir','out')
c.TickDirection = 'out'; c.Label.String = 'cells/250';
ylim([13 22])

nexttile
VS = mean(cat(4,DCS.Result),4);
imagesc(VS')
c = colorbar;
clim([0 20])
title('OFC-CPr')
set(gca,'XColor', 'none','YColor','none','box','off','Colormap', cmap3,'TickDir','out')
c.TickDirection = 'out'; c.Label.String = 'cells/250\mum';
ylim([13 22])

nexttile
axis off
nexttile
axis off
nexttile
axis off
nexttile
axis off
nexttile
axis off
nexttile
axis off

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



