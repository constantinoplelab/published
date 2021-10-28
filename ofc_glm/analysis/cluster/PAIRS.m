function [pairs,p_pairs,pairs_ref,meanvals_dat,meanvals] = ...
        PAIRS(fet_proc,kval,nset,saveops)
%calculates projection angle of response similarity statistic
%
%INPUTS:
%   fet_proc: (nn x dim) features space for nn neurons, dim features.
%       should processed features with pca to have zero mean and diag cov.
%   kval: number of nearest neighbors in averaging
%   nset: number of reference distributions to build
%   saveops: ops struct for how handle saving/loading reference data
%       saveref: boolean. if true, save reference data
%       loadref: boolean. if true, will load refence data instead of calc
%       savename: file name for where to save reference data
%       loadname: file name from where to loaad reference data
%
%OUTPUTS:
%   pairs: PAIRS statistic value
%   p_pairs: p value of PAIRS stat, based on norm. dist of ref data medians
%   pairs_ref: PAIRS statistic for grand median of ref data and each ref.
%       dataset
%   meanvals_dat: (nn x 1) value of average angle to nearest neighbors
%   meanvals: (nn x nset) values of anverage angles for each ref. dataset


%function for angle
thetafun = @(u,v) real(acos(max(min(dot(u,v)/(norm(u)*norm(v)),1),-1)));


%% whiten features
nn = size(fet_proc,1);
dim = size(fet_proc,2);
Sigma = cov(fet_proc);
W = real(inv(Sigma).^(0.5)); %covariance was diagonal from projection on PCA. square root is elementwise

fet_whitened = fet_proc*W;

%% build/load reference data

if ~saveops.loadref
    disp('building reference data')
    fet_ref = zeros(nset,nn,dim);
    for j = 1:nset
        fet_ref(j,:,:) = mvnrnd(zeros(dim,1),eye(dim),nn);
    end
  
else 
    disp('loading previously calc. reference data')
    f = load(saveops.loadname);
    fet_ref = f.fet_ref;
    meanvals = f.meanvals; %average angle per sample
end
    

if saveops.saveref
    disp('saving ref data')
    save(saveops.savename,'fet_ref')
end


%% LONG CALC: calculate average theta vals per neighbor


if ~saveops.loadref
    disp('calculating angles in ref data. long calc...')
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

end


%when done delete the reference dist
if saveops.saveref
    save(saveops.savename,'meanvals','-append');
    clear fet_ref f
end


%%  calculate the mean angle for true data

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

%% calculate PAIRS

%true PAIRS for data
thetamedian_ref_all = median(meanvals);
pairs = (thetamedian_ref_all - thetamedian_dat)/thetamedian_ref_all;

%PAIRS for refernce data
pairs_ref = zeros(nset,1);
for j = 1:nset
    dat_ref = squeeze(meanvals(j,:));
    thetamedian_ref = median(dat_ref(:));

    pairs_ref(j) = (thetamedian_ref_all - thetamedian_ref)/thetamedian_ref_all;
end


%either 1.95 STD, or numerical 95% two-sided CI
mu_ref = mean(pairs_ref);
std_ref = std(pairs_ref);

p_pairs = 2*(normcdf(pairs,mu_ref,std_ref,'upper'));


