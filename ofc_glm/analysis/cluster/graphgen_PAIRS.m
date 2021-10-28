%study to look at PAIRS statistic, via 2014 nat neuro paper

addpath(genpath('~/projects/ofc'))
%datadir =  '/Users/dhocker/projects/ofc/results/clustering/PSTH1/Zbs90_ns100pca95/';
datadir =  '/Users/dhocker/projects/ofc/results/clustering/singleConditional2/Zbs90_ns100pca95/';
f = load(strcat(datadir,'features.mat'));

refdat_file = strcat(datadir,'pairs_ref.mat'); %where raw distributions of reference data are stored
meanvals_ref_fname = strcat(datadir,'pairs_meanvals_ref.mat'); %where mean values of angles for ref are stored

%kval = 3; % for PSTH clustering
kval = 8; %for conditional clustering

%function for angle
thetafun = @(u,v) real(acos(max(min(dot(u,v)/(norm(u)*norm(v)),1),-1)));

savedat = false;
loaddat = true;

nset = 10000;
%nset = 100; %for checking k value

%% examine featuers. they are zero-meaned, but not whitened

fet_proc = f.fet_proc(:,:,1);

rng(101)

nn = size(fet_proc,1);
dim = size(fet_proc,2);



Sigma = cov(fet_proc);
W = real(inv(Sigma).^(0.5)); %covariance was diagonal from projection on PCA. square root is elementwise

fet_whitened = fet_proc*W;


%% build reference data


fet_ref = zeros(nset,nn,dim);
for j = 1:nset
    fet_ref(j,:,:) = mvnrnd(zeros(dim,1),eye(dim),nn);
end

savedat = false;
if savedat    
    save(refdat_file,'fet_ref');
end




%% FIND BEST K calcualte all pairwise angles. THIS IS SLOW. do here for nset=100 ONLY

thetavals = zeros(nset,nn,nn);

for j = 1:nset
    disp(j)
    for k = 1:nn
        for m = k:nn
            thetavals(j,k,m) = thetafun(squeeze(fet_ref(j,k,:)),squeeze(fet_ref(j,m,:)));
            thetavals(j,m,k) = thetavals(j,k,m);
        end
    end
end


%% FIND BEST K.  for a given k, calcualte the average angle between neighbors

kvec = 2:20;
setind = 1; %just look at one dataset for now


median_theta_dist = zeros(nset,numel(kvec)); %median k-neighbor aveage from dis.
meanvals = zeros(nset,numel(kvec),nn); %the k-neighbor average for each neuron, given k

for m = 1:nset
    disp(m)
    dat_k = squeeze(thetavals(m,:,:));   
    for j = 1:numel(kvec)        
        for k = 1:nn
            [~,I] = sort(dat_k(k,:),'ascend');
            meanvals(m,j,k) = mean(dat_k(k,I(2:kvec(j)+1)));
        end  

        median_theta_dist(m,j) = median(meanvals(j,:));
    end
end

%% FIND BEST K visualize the distribution for different k and the median
figure(1)
subplot(1,3,1)
bins = linspace(0,pi,100);

ind = 1;
dat_h = squeeze(meanvals(:,ind,:));
histogram(dat_h(:),bins,'Normalization','probability')
vline(median(dat_h(:)),'r')
vline(pi/4,'k')
title(strcat('k = ',num2str(kvec(ind))))
xlabel('angle (rad)')
xlim([0.5,1.5])
ylabel('p')
%legend('dist','median','pi/4')

subplot(1,3,2)
ind = 6;
dat_h = squeeze(meanvals(:,ind,:));
histogram(dat_h(:),bins,'Normalization','probability')
vline(median(dat_h(:)),'r')
vline(pi/4,'k')
title(strcat('k = ',num2str(kvec(ind))))
xlabel('angle (rad)')
xlim([0.5,1.5])
ylabel('p')
%legend('dist','median','pi/4')

subplot(1,3,3)
ind = 7;
dat_h = squeeze(meanvals(:,ind,:));
histogram(dat_h(:),bins,'Normalization','probability')
vline(median(dat_h(:)),'r')
vline(pi/4,'k')
title(strcat('k = ',num2str(kvec(ind))))
xlabel('angle (rad)')
xlim([0.5,1.5])
ylabel('p')
%legend('dist','median','pi/4')

%% LONG CALC: after best k is chosen, recalculate mean values and discard thetavals

%load data
if loaddat
    g = load(refdat_file);
    fet_ref = g.fet_ref;
end

kval = 3;
meanvals = zeros(nset,nn); %the k-neighbor average for each neuron, given k

for m = 1:nset
    disp(m)
    %calculate thetavals from the reference distribution for that point
    thetavals_m = zeros(nn);
    for k = 1:nn
        for j = k:nn
            thetavals_m(k,j) = thetafun(squeeze(fet_ref(m,k,:)),squeeze(fet_ref(m,j,:)));
            thetavals_m(j,k) = thetavals_m(k,j);
        end
    end
        
    for k = 1:nn
        [~,I] = sort(thetavals_m(k,:),'ascend');
        meanvals(m,k) = mean(thetavals_m(k,I(2:kval)));
    end  
end



%when done delete the reference dist
if savedat 
    save(meanvals_ref_fname,'meanvals');
    clear fet_ref, g
end



%% REPLOT, STARt HERE with a chosen K, calculate the PAIRS statistic

%this was moved to graphgen_PAIRS

%kval = 3;
%kval = 8

thetavals_dat = zeros(nn);

for k = 1:nn
    for m = k:nn
        thetavals_dat(k,m) = thetafun(squeeze(fet_whitened(k,:)),squeeze(fet_whitened(m,:)));
        thetavals_dat(m,k) = thetavals_dat(k,m);
    end
end

meanvals_dat = zeros(nn,1); %the k-neighbor average for each neuron, given k


dat_k = thetavals_dat;  
for k = 1:nn
    [~,I] = sort(dat_k(k,:),'ascend');
    meanvals_dat(k) = mean(dat_k(k,I(2:kval)));
end  


thetamedian_dat = median(meanvals_dat);

%% calculate the distribution of this test statistic for the datasets

if loaddat
    g = load(meanvals_ref_fname);
    meanvals = g.meanvals;
end

%true PAIRS for data
thetamedian_ref_all = median(meanvals);
pairs_all = (thetamedian_ref_all - thetamedian_dat)/thetamedian_ref_all;

%PAIRS for refernce data
pairs_ref = zeros(nset,1);
for j = 1:nset
    dat_ref = squeeze(meanvals(j,:));
    thetamedian_ref = median(dat_ref(:));

    pairs_ref(j) = (thetamedian_ref_all - thetamedian_ref)/thetamedian_ref_all;
end


figure(2)
clf
histogram(pairs_ref,100,'Normalization','probability')
CI_95 = 1.95*std(pairs_ref);
vline(mean(pairs_ref)+CI_95,'k')
vline(mean(pairs_ref)-CI_95,'k')
vline(pairs_all,'r')
xlabel('PAIRS')
title('Reference data PAIRS comparison to data PAIRS statistic (red)')

%calculate a confidence interval around this distribution

%either 1.95 STD, or numerical 95% two-sided CI
mu_ref = mean(pairs_ref);
std_ref = std(pairs_ref);

%p = 2*normcdf(pairs_all,mu_ref,std_ref)
p = 2*(normcdf(pairs_all,mu_ref,std_ref,'upper'))




%% make a simplified plot showing angle distribution for real data vs grand reference data

figure(3)
clf
hold on

bins = linspace(0,pi/3,100);
histogram(meanvals(:),bins,'Normalization','probability')
hold on
histogram(meanvals_dat,bins,'Normalization','probability')
title('Distribution of angles for population, PSTH clustering')
xlabel('angle (rad)')
xlim([0,pi/3])
ylabel('p')
legend('reference data','true data')

vline(median(meanvals_dat),'k')
vline(median(meanvals(:)),'k')


[h,p] = kstest2( meanvals(:) , meanvals_dat )


