%recreates clustering results of Figure 1, on PSTH-based feature space

datadir = '/Users/dhocker/projects/ofc/data/published/';


%% load cluster labels and aggregated data structure

disp('loading aggregated data and cluster labels');

clustfilename = strcat(datadir,'clustlabs_psth.mat');
clustfile = load(clustfilename);
clustlabs = double(clustfile.labels)'-1;


clustfilename_cond = strcat(datadir,'clustlabs_cond.mat');
clustfile_cond = load(clustfilename_cond);
clustlabs_cond = double(clustfile_cond.labels)'-1;
%swap 0 and 1 for better interpretability across feature spaces
clustlabs_cond(clustlabs_cond==0) = 10;
clustlabs_cond(clustlabs_cond==1) = 0;
clustlabs_cond(clustlabs_cond==10) = 1;

f = load(strcat(datadir,'concatdata_ofc_noOptOut.mat'));
A = f.A;
%clear f


%% load PSTHs into a single data structure

disp('calculating PSTHs')

nn = 659; %number neurons
nt = 121; %number datapoints


psth = zeros(nn,nt);

for j = 1:nn
    psth(j,:) = nanmean(A{j}.hmat_start);
    
end



%% plot 

disp('plotting')
useZscore = true;

[fig_psth,fig_psthheatmat] = ...
    graphgen_psth_cluster(psth,clustlabs,clustfilename,useZscore);

set(fig_psth,'Position',[29   378   560   420]);
set(fig_psthheatmat,'Position',[595   381   560   420]);


%% plot conditional 

disp('plotting')
useZscore = true;

[fig_cond,fig_condheatmat] = ...
    graphgen_psth_cluster(psth,clustlabs_cond,clustfilename_cond,...
    useZscore);

set(fig_cond,'Position',[29   978   560   420]);
set(fig_condheatmat,'Position',[595   981   560   420]);
