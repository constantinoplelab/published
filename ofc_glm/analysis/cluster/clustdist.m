function [Dmat,mumat,S] = clustdist(fet,labs)
%calculate mahlanobis distances of points in feature space to cluster distributions


nclust = numel(unique(labs));
[nn,nx] = size(fet); %number samples, number of features

S = zeros(nclust,nx,nx); %cluster covariance
mumat = zeros(nclust,nx); %cluster means
Dmat = zeros(nn,nclust); %distances to each cluster distribution


for j = 1:nclust
    Sj = cov(fet(labs==j-1,:));
    mj = mean(fet(labs==j-1,:))';
    Sjinv = inv(Sj);
    
    for k = 1:nn
        Dmat(k,j) = sqrt((fet(k,:)'-mj)'*Sjinv*(fet(k,:)'-mj));
    end
    
    S(j,:,:) = Sj;
    mumat(j,:) = mj;
end

%% plot

distmat = zeros(nclust);
for j = 1:nclust
    for k = 1:nclust
        distmat(j,k) = nanmean(Dmat(labs==j-1,k));
    end
end
        
figure()
clf
imagesc(distmat)
xlabel('cluster distribution, B')
ylabel('cluster A')
title('average distribution distance E_A[D(x_A,B)]')
xticks(1:nclust)
yticks(1:nclust)
colorbar()
%caxis([0,6])
