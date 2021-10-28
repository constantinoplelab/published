function [dp,dpshuff,dp_raw] = dprime(hmat,mask1,mask2,varargin)
%dprime: discriminability in firing rate between two conditions
%
%INPUTS:
%   hmat: (ns, nt) heatmat for a single neuron.ns trials. nt time points
%   mask1: (ns,1) boolean for condition 1. true if trial satistfies cond.
%   mask2: similarly to mask1. Typically mask1 and mask2 should not share
%           true values
%
%OUTPUTS:
%   dp: (nt,1) shuffle mean-subtracted, sig.-corrected discriminability 
%   dpshuff: (nt,nshuff) samples of shuffled dprime values. raw data used 
%       to calc  shufflemean and confidence interval
%   dp_raw(nt,1) original dprime values, before mean subtraction and sig.
%       correciton

%% variable arguments defined here, so we don't have to dig for them
if isempty(varargin)
    ci = 0.95; %confidence interval
    useZscore = true; %z-score responses before taking dprime? don't think is matters.
else
    switch numel(varargin)
        case 1
            ci = varargin{1}; %confidence interval
        case 2
            ci = varargin{1};
            useZscore = varargin{1};
    end
    
end

%% calculate

[~,nt] = size(hmat); 
rng(101) %reproducible shuffles


mu0 = nanmean(nanmean(hmat));
std0 = nanstd(nanmean(hmat));
if useZscore
    dat1 = (hmat(mask1,:)-mu0)/std0;
    dat2 = (hmat(mask2,:)-mu0)/std0;
else
    dat1 = hmat(mask1,:);
    dat2 = hmat(mask2,:);
end


%shuffle setup: balance conditions by upsampling (sampling with replacement)
nsfmax = max([sum(mask1),sum(mask2)]); %condition w largest # trials
trials1 = find(mask1);
trials2 = find(mask2);
trials_shuff = sort([randsample(trials1,nsfmax,true);randsample(trials2,nsfmax,true);]);
nsf = nsfmax*2; %total trials in data after balancing conditions

%create shuffled trial numbers for each sample of the shuffle dist.
nshuff = 1000;
shuffmat = zeros(nshuff, nsf);
for j = 1:nshuff
    shuffmat(j,:) = randperm(nsf);
end

dpshuff = zeros(nshuff,nt); %TODO: should be dpshuff

%make the raw PSTH data for the shuffle, and equivalent condition masks
if useZscore
    hmatbal = (hmat(trials_shuff,:)-mu0)./std0;
else
    hmatbal = hmat(trials_shuff,:);
end
mask1bal = mask1(trials_shuff); 
mask2bal = mask2(trials_shuff);

%d', true data
dp_raw = abs(nanmean(dat1)-nanmean(dat2))./sqrt(0.5*(nanvar(dat1)+nanvar(dat2)));

%d' for shuffle
for m = 1:nshuff

    %shuffle the masks
    mask1shuff = mask1bal(shuffmat(m,:));
    mask2shuff = mask2bal(shuffmat(m,:));

    %get scrambled data
    dat1 = hmatbal(mask1shuff,:);
    dat2 = hmatbal(mask2shuff,:);

    %d', shuffle
    dpshuff(m,:) = abs(nanmean(dat1)-nanmean(dat2))./sqrt(0.5*(nanvar(dat1)+nanvar(dat2)));     

end

% signficance correct
shuffmean = squeeze(nanmean(dpshuff));
%shuffstd = squeeze(nanstd(dpshuff));

%find numerical one-sided CI
shuffsig = zeros(size(dp_raw));%cutoff value of the confidence interval
for m = 1:nt
    [h,x] = histcounts(dpshuff(:,m),20,'Normalization','cdf');
    sig_idx = find(h > ci);
    if ~isempty(sig_idx)
        shuffsig(m) = x(sig_idx(1)); %first instance
    end
end

%subtract mean of shuffle and mask out whatever is below confidence
%interval
dp = dp_raw-shuffmean;
dp(dp_raw-shuffsig < 0) = 0;
%dtmat_sig(dtmat-shuffsig < 0) = nan;


    

    


