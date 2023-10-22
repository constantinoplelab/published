function [shuffleTauDiff, observedTauDiff, pval] =...
    myTauShuffleTest(first_q, last_q, myBetas, Nshuffle)
Betas = myBetas(first_q | last_q,:);
nrats = size(Betas, 1);

shuffleTau = nan(Nshuffle, 2);

for ii = 1:Nshuffle
    fprintf('shuffle %d out of %d\n', ii, Nshuffle)
    shuffleLatency = randperm(nrats)';
    
    group1 = shuffleLatency(1:floor(nrats)/2);
    Betas1 = mean(Betas(group1,1:end-1));
    [~, bestFit1] = fit_exp_decay(1:10, 10*Betas1, 10,...
        [-3 0], [1 50]);
    shuffleTau(ii,1) = bestFit1(2);
    
    group2 = shuffleLatency(floor(nrats)/2+1:end);
    Betas2 = mean(Betas(group2,1:end-1));
    [~, bestFit2] = fit_exp_decay(1:10, 10*Betas2, 10,...
        [-3 0], [1 50]);
    shuffleTau(ii,2) = bestFit2(2);
end

latencyBetasObsv1 = mean(myBetas(first_q,1:end-1));
[~, bestFitObsv1] = fit_exp_decay(1:10, latencyBetasObsv1, 10,...
    [-10 0], [1 10]);

latencyBetasObsv2 = mean(myBetas(last_q,1:end-1));
[~, bestFitObsv2] = fit_exp_decay(1:10, latencyBetasObsv2, 10,...
    [-10 0], [1 10]);

shuffleTauDiff = shuffleTau(:,1) - shuffleTau(:,2);
observedTauDiff = bestFitObsv1(2) - bestFitObsv2(2);

pval = sum(abs(shuffleTauDiff) > abs(observedTauDiff))/(2*Nshuffle);
end
