function plot_fig6h(datadir)
% Plot average cross-validated classification accuracy of a logistic
% regression decoder. N = 18 sessions from 8 rats.
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

fpath = fullfile(datadir, 'data-published/NpxData/npx_DMS.mat');
load(fpath, 'S', 'SU');
[cvPerf, ~] = decodeSpeed(SU, S);

figure
set(gcf, units='inches', renderer='painter', ...
    position=[9,5,1.6,1.9])
tiledlayout('flow');
nexttile;
hold on
jitx = rand(size(cvPerf,1),1)*0.1;
jitx(1:end/2) = -jitx(1:end/2);
plot(1+jitx, cvPerf(:,1), 'ko')
errorbar(1, mean(cvPerf(:,1), 'omitnan'), std(cvPerf(:,1), 'omitnan'), ...
    'k_', capsize=0, linewidth=1.5)
yline(50, 'k--')
xlim([0.5 1.5])
ylabel('Cross-validated accuracy (%)')
xticks([])
set(gca, tickdir='out', xtick=1, fontsize=8)

end

function [cvPerf, ridgeWeights] = decodeSpeed(SU, S)
% calculate cross-validated logistic regression linear classifier error
% of predicting fast vs slow trials from firing rate at Reward Port On.
% Contralateral trials only.


twin = [0 0.5];
cvErr = nan(length(SU),1);
ridgeWeights = nan(length(SU),1);

for sess=1:length(S)
    if length(SU{sess}) < 5 % exclude sessions with fewer than 5 cells
        continue
    else
        ratName = S{sess}.RatName;

        % get contralateral trials
        if sum(strcmpi(ratName, {'j081', 'j076'}))>0
            recHemi = 'R';
        else
            recHemi = 'L';
        end
        isContra = cellfun(@(x) ~strcmpi(x, recHemi), S{sess}.RewardedSide)';

        % For Reward Port On:
        % get reaction time 
        RT = S{sess}.slrt;
        % filter outliers
        thresh = [];
        thres(1) = prctile(RT, 1); % lower threshold
        thresh(2) = prctile(RT, 99); % upper threshold
        RT(RT<=thresh(1)) = nan;
        RT(RT>=thresh(2)) = nan;
        % nan ipsi trials
        RT(~isContra) = nan;

        % bin reaction times
        quantiles = [0 1/4 2/4 3/4 1];
        bins = quantile(log2(RT), quantiles);
        
        % get fast vs slow trials
        isFast = log2(RT)>=bins(1) & log2(RT)<=bins(2);
        isSlow = log2(RT)>=bins(end-1) & log2(RT)<=bins(end);
        fprintf('Slow: %i trials, Fast: %i trials\n', ...
            sum(isSlow), sum(isFast))
                
        % create true/false vector 
        y = find(isFast | isSlow);
        y(ismember(y,find(isFast))) = 1;
        y(ismember(y,find(isSlow))) = 0;
            
        % construct trial x neuron x T tensor
        X = [];

        T = SU{sess}{1}.xvec.SON;
        [~,t1] = min(abs(T-twin(1)));
        [~,t2] = min(abs(T-twin(2)));
        
        for n=1:length(SU{sess})
            these = isFast | isSlow;
            X(:,n,:) = SU{sess}{n}.hmat.SON(these,t1:t2);
        end

        % flatten
        X = reshape(X, sum(these), []);

    
        % fit model
        Lambda = logspace(-4,-0.2,11);
        mdl = fitclinear(X', y, ...
            'ObservationsIn', 'columns', ...
            'CrossVal', 'on', 'Learner', 'logistic', ...
            'Regularization', 'ridge', 'Lambda', Lambda);
        [cvErr(sess,1), ridgeWeights(sess,1)] = min(kfoldLoss(mdl));

    end
end

%
rmv = isnan(cvErr(:,1));
cvErr(rmv,:) = [];
ridgeWeights(rmv,:) = [];
cvPerf = 100*(1-cvErr);


end




