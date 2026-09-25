% figure generating scripts for probability simplex
% and sample neuron PSTHS, conditioned on most likely or 2nd most likely
% block


%% add path and load DS projection data

%location of analysis codebase
codepath = '~/david/maggie';

% add code paths
addpath(genpath(codepath))
addpath(genpath(codepath + "david"))
addpath(genpath(codepath + "maggie"))

%for saving pds correctly so illustrator can modify them
set(0, 'DefaultFigureRenderer', 'painters');

% DS projection neuron data location
%datadir = "/Users/dhocker/projects/dynamics/data/maggie/";
datadir = "Z:\david\maggie\Most_Likely_block_analysis\";
fname = datadir + "DS_projection_neurons_non-Stimulated.mat";
b = load(fname);

%% run bayes model for each session and get MLB, MLB2
% for every cell (only 60 now, and 23 sessions, so not prohibitive), do bayes calc
% bayes params for all sessions
kappa_mi = 1; %.5;  
kappa_hi = 1.5; %1.2; %0.7; 
kappa_lo = 0.5;%0.8; %.3; 
D = 0.5; 
lambda =1;
noise = 80;
params = [kappa_mi, kappa_hi, kappa_lo, D, lambda];
nsess  = numel(b.all_S_DS);

Beliefs_all = [];
wait_time_all = [];
Beliefs_bysess = {};
for j = 1:nsess
    disp(j)
    %get session data
    S_j = b.all_S_DS{j};

    % this struct is almost the correct format for an A struct, but needs
    % a little renaming to run. 
    S_j.ntrials = [numel(S_j.Block)]; % only 1 session in this struct
    S_j.reward = S_j.RewardAmount;
    S_j.prob_catch = S_j.ProbCatch;

    [~, wait_time, ~, Belief, ~, ~] = GenerateSynthData_Bayes(params, S_j,'logn', true,noise);
    Beliefs_bysess{j} = Belief;
    Beliefs_all = [Beliefs_all, Belief];
    wait_time_all = [wait_time_all, wait_time'];

end

%remove nans
Beliefs_m = Beliefs_all(1,:);
Beliefs_h = Beliefs_all(2,:);
Beliefs_l = Beliefs_all(3,:);
Beliefs_m(isnan(wait_time_all)) = [];
Beliefs_h(isnan(wait_time_all)) = [];
Beliefs_l(isnan(wait_time_all)) = [];

Beliefs_scrubbed = [Beliefs_m;Beliefs_h; Beliefs_l];
    
%calculate mostlikely block (MLB) and 2nd most likely block (MLB2)
[B,I] = sort(Beliefs_scrubbed,1,'descend');
mlb = I(1,:);
mlb2 = I(2,:);

%% plot simplex
figure(58)
clf
hold on

grid
xlabel('p mix')
ylabel('p high')
zlabel('p low')
title('belief distribution')
set(gca,'fontsize',15)

 
x = [ 1, 0.5, 1/3, 0.5];
y = [ 0, 0.5, 1/3, 0];
z = [ 0, 0,   1/3, 0.5];
fill3(x,y,z,'k','facealpha',0.2)

%high
x = [ 0, 0.5, 1/3, 0.0];
y = [ 1, 0.5, 1/3, 0.5];
z = [ 0, 0,   1/3, 0.5];
fill3(x,y,z,'r','facealpha',0.2)

%low
x = [ 0, 0.0, 1/3, 0.5];
y = [ 0, 0.5, 1/3, 0.0];
z = [ 1, 0.5,   1/3, 0.5];
fill3(x,y,z,'b','facealpha',0.2)
view(69,25)

% plot beliefs, color coded by MLB2. give noise jitter to location for vis.
ntrial = size(Beliefs_scrubbed,2);
noise = 0.005*randn(3,ntrial);
scatter3(Beliefs_scrubbed(1,mlb2==1)+noise(1,mlb2==1), Beliefs_scrubbed(2,mlb2==1)+noise(2,mlb2==1), Beliefs_scrubbed(3,mlb2==1)+noise(3,mlb2==1),30,'o','markerfacecolor','k','markeredgecolor','k','markerfacealpha',0.2,'markeredgealpha',0.2)
scatter3(Beliefs_scrubbed(1,mlb2==2)+noise(1,mlb2==2), Beliefs_scrubbed(2,mlb2==2)+noise(2,mlb2==2), Beliefs_scrubbed(3,mlb2==2)+noise(3,mlb2==2),30,'o','markerfacecolor','r','markeredgecolor','r','markerfacealpha',0.2,'markeredgealpha',0.2)
scatter3(Beliefs_scrubbed(1,mlb2==3)+noise(1,mlb2==3), Beliefs_scrubbed(2,mlb2==3)+noise(2,mlb2==3), Beliefs_scrubbed(3,mlb2==3)+noise(3,mlb2==3),30,'o','markerfacecolor','b','markeredgecolor','b','markerfacealpha',0.2,'markeredgealpha',0.2)

ylim([-0.0,1])
xlim([-0.0,1])
zlim([-0.0,1])

%% plot MLB or MLB2 in select neurons

usemlb = false; % if true, most likely block. if false, 2nd mostl likely

%choose a cell. the list below are thigns i found from regression analysis
%that had strong regression weights in certain areas. not an exhuastie
%list, but good place to start -------

%cell_idx = 42; % mlb
%cell_idx = 51; %mlb 
%cell_idx = 2; %mlb2
%cell_idx = 40; %mlb2

% cells good at 1s before reward.
%cell_idx = 47;
%cell_idx = 16;
%cell_idx = 45;
%cell_idx = 50;
%cell_idx = 2;
%cell_idx = 9;
%cell_idx = 35;
%cell_idx = 52;

% for MLB: 1,4, 31, 45(Big), 49(big), 60(big), 62(big), 79, 94(big), 101,
% 114
%cell_idx = 3;
%cell_idx = 34;
%cell_idx = 48;
%cell_idx = 52;
%cell_idx = 2;
%cell_idx = 4;
%cell_idx = 21;
%cell_idx = 40;
cell_idx = 47;
%cell_idx = 60;

%cell_idx = 11;
%cell_idx = 52;

%find its session
sess_idx = b.all_index_DS{cell_idx,4};
xvec = b.all_SU_DLS_200{1}.xvec.Rew;

[B,I] = sort(Beliefs_bysess{sess_idx},1,'descend');
mlb_neuron = I(1,:);
mlb2_neuron = I(2,:);

%choose mlb or mlb2

if usemlb
    cond2use = mlb_neuron;
else
    cond2use= mlb2_neuron;
end

% get data, calculate mean and sem for condition
dat = b.all_SU_DLS_200{cell_idx}.hmat.Rew;
hits = b.all_S_DS{sess_idx}.hits;

figure(168)
clf
hold on

mask = hits'==1 & cond2use ==2;
nm = sum(mask);
dmean = mean(dat(mask,:),1,'omitnan');
dsem = std(dat(mask,:),[],1,'omitnan')/sqrt(nm);
shadedErrorBar(xvec,dmean,dsem,'lineprops',{'color',[0.8,0,0.0]} )

mask = hits'==1 & cond2use ==1;
nm = sum(mask);
dmean = mean(dat(mask,:),1,'omitnan');
dsem = std(dat(mask,:),[],1,'omitnan')/sqrt(nm);
shadedErrorBar(xvec,dmean,dsem,'lineprops',{'color',[0.5,0,0.8]} )

mask = hits'==1 & cond2use ==3;
nm = sum(mask);
dmean = mean(dat(mask,:),1,'omitnan')
dsem = std(dat(mask,:),[],1,'omitnan')/sqrt(nm)
shadedErrorBar(xvec,dmean,dsem,'lineprops',{'color',[0.0,0,0.8]} )

vline(0,'k')
xlabel('Time to reward (s)')
xlim([-1 2.5])
ylabel('firing rate (Hz)')

if usemlb
    title(strcat('neuron ',num2str(cell_idx),'most likely block encoding'))
    legend('mlb = high','mlb = mixed','mlb = low')
else

    title(strcat('neuron ',num2str(cell_idx),', 2nd most likely block encoding'))
    legend('mlb2 = high','mlb2 = mixed','mlb2 = low')
end

set(gca,'fontsize',15,'Box','off','TickDir','out')
xticks(-1:1:3)



mask = hits'==1 & cond2use ==2;
nm = sum(mask);
d1 = dat(mask,:);
mask = hits'==1 & cond2use ==1;
nm = sum(mask);
d2 = dat(mask,:);
mask = hits'==1 & cond2use ==3;
nm = sum(mask);
d3 = dat(mask,:);

d4 = ones(1,61)*max([d1;d2;d3],[],'all')*2;

load('Z:\david\maggie\CBar1.mat')
if usemlb==true
figure(169)
imagesc(xvec,linspace(1,length([d1;d4;d2;d4;d3])), [d1;d4;d2;d4;d3])
cb = colorbar;
colormap(CBar1)
cb.Ticks = 0:10:30;
clim([0 40])
set(gca,'yTick',(0:20:400),'Box','off','TickDir','out')
xlim([-1 3])
xticks(-1:1:2.5)
xlabel('Time from reward')
xline(0,'--w')
cb.Label.String = 'Firing rates (Hz)';
ylabel('Trials')
else
%     figure(169)
% imagesc(xvec,linspace(1,length([d1;d4;d2;d4;d3])), [d1;d4;d2;d4;d3])
% cb = colorbar;
% colormap(CBar1)
% cb.Ticks = 0:10:30;
% clim([0 40])
% set(gca,'yTick',(0:20:400),'Box','off','TickDir','out')
% xlim([-1 3])
% xticks(-1:1:2.5)
% xlabel('Time from reward')
% xline(0,'--w')
% cb.Label.String = 'Firing rates (Hz)';
% ylabel('Trials')

  figure(169)
imagesc(xvec,linspace(1,length([d1;d4;d2;d4;d3])), [d1;d4;d2;d4;d3])
cb = colorbar;
colormap(CBar1)
cb.Ticks = 0:10:30;
clim([0 30])
set(gca,'yTick',(0:20:400),'Box','off','TickDir','out')
xlim([-1 3])
xticks(-1:1:2.5)
xlabel('Time from reward')
xline(0,'--w')
cb.Label.String = 'Firing rates (Hz)';
ylabel('Trials')
end


bins = 0.2;
raster = makeHeatmat(all_SU_DS(cell_idx), all_S_DS{all_index_DS{cell_idx,4}}, all_S_DS{all_index_DS{cell_idx,4}}.behEvents, [-4 8],bins,1);
dat2 = raster{1}.raster.Rew;

figure

subplot(3,1,1)
mask = hits'==1 & cond2use ==2;
dat2a = dat2(mask);
for i = 1:length(dat2a)
    if ~isempty(dat2a{i})
y = [];
y = ones(length(dat2a{i}))*i;
plot(dat2a{i},y,'|k')
hold on
    end
end


subplot(3,1,2)
mask = hits'==1 & cond2use ==1;
dat2a = dat2(mask);
for i = 1:length(dat2a)
    if ~isempty(dat2a{i})
y = [];
y = ones(length(dat2a{i}))*i;
plot(dat2a{i},y,'|k')
hold on
    end
end


subplot(3,1,3)
mask = hits'==1 & cond2use ==3;
dat2a = dat2(mask);
for i = 1:length(dat2a)
    if ~isempty(dat2a{i})
y = [];
y = ones(length(dat2a{i}))*i;
plot(dat2a{i},y,'|k')
hold on
    end
end

labelname = {'High','Mix','Low'};
for s = 1:3
    subplot(3,1,s)
    xlim([-1 2.5]); xticks(-1:1:3)
xline(0,'--k')
box off
set(gca,'TickDir','out')
yticks([])
ylabel(labelname{s})
end

if usemlb
    sgtitle('Most likely block');
else
    sgtitle('Second most likely block');
end




