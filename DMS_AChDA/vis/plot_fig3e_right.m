function plot_fig3e_right(datadir)
% Plot dopamine AUC at the offer cue (0-0.5 s) vs -delta trial initiation
% time, averaged across rats. N = 10 rats.
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

fname = fullfile(datadir, 'data-published/deltaITI_vs_AUC.mat');
if exist(fname, 'file')
    load(fname, 'auc', 'ratList')
else
    % load all DA rats
    ratList = loadRats(datadir, 'da');    
    auc = cell(length(ratList), 1);
    bins = auc;
    for rr=1:length(ratList)
        ratname = ratList{rr};
        fprintf('%i of %i: %s\n', rr, length(ratList), ratname)
        [auc{rr}, bins{rr}] = deltaITI_vs_AUC(datadir, ratname);
    end
    auc = cell2mat(auc);
    save(fname, 'auc', 'bins', 'ratList');
end

nbin = size(auc,2);

% normalization
auc_norm = (auc - mean(auc,2)) ./ std(auc,0,2);

% stats
nRats = length(ratList);

ratID = repmat((1:nRats)', nbin, 1);
binNum = repelem((1:nbin)', nRats);
valVec = auc_norm(:);

tbl = table(categorical(ratID), binNum, valVec, ...
    'VariableNames', {'ratID','binNum','valVec'});

% Fit random-intercept model
lme = fitlme(tbl, 'valVec ~ binNum + (1|ratID)');

disp(lme)

% One-sided p-value for slope > 0
beta = fixedEffects(lme);
se   = lme.Coefficients.SE(2); % slope SE
df   = lme.DFE;
tval = beta(2)/se;
p_one_sided = 1 - tcdf(tval, df);

fprintf('Slope = %.4f, one-sided p = %.4g\n', beta(2), p_one_sided);

figure; 
set(gcf, units='inches', renderer='painters',...
    position=[7,3,2,2])
plotPretty(1:nbin, auc_norm, 'k');
xlim([0.5,nbin+0.5])
xticks(1:nbin)
ylim([-1.5 1])
yticks([-1:1:1])
yline(0, '--', color='#a0a0a0')
xline(3.5, '--', color='#a0a0a0')

% fit a line
lm = fitlm(1:nbin, mean(auc_norm, 1));
b = lm.Coefficients.Estimate;
hold on
xx = linspace(0, 9, 10);
plot(xx, xx*b(2)+b(1), '-', color='#6D6E71', linewidth=1)

xlabel('-\Delta ITI bin')
ylabel('DA AUC')
set(gca, fontsize=10)
subtitle(sprintf('N = %i rats', length(ratList)))

end


function [auc_mean, edges] = deltaITI_vs_AUC(datadir, rat)
% CPIn AUC on trial N vs -(ITI(N+1)-ITI(N)) = ITI(N)- ITI(N+1)
% big RPE at CPIn should cause a big change in ITI 

% load pstruct and bstruct
basefolder = fullfile(datadir, 'data-published', 'PhotometryData');
% DMS DA rats are either in GRAB_DA_DMS or GRAB_rDAgACh_DMS folders
fname = strcat(rat, '*DA_ch*_DMS.mat');
pfile = dir(fullfile(basefolder, 'GRAB_DA_DMS', fname));
bfile = dir(fullfile(basefolder, 'GRAB_DA_DMS', ...
    strcat(rat, '_DA_bstruct.mat')));
if isempty(pfile)
    pfile = dir(fullfile(basefolder, 'GRAB_rDAgACh_DMS', fname));
    bfile = dir(fullfile(basefolder, 'GRAB_rDAgACh_DMS', ...
        strcat(rat, '_rDAgACh_bstruct.mat')));
end
load(fullfile(pfile.folder, pfile.name), 'pstruct')
pdata = pstruct.(rat).CPIn;
load(fullfile(bfile.folder, bfile.name), 'bstruct')
bdata = bstruct.(rat).CPIn;

% nan ITIs on first trial
bdata.ITI(bdata.TrialNumber==1) = nan;

% determine iti cutoffs using violation and post-violation trials
vios = bdata.TrialType==2;
postvios = [0;vios(1:end-1)];

cutoff = [prctile(bdata.ITI(vios|postvios),10), ...
    prctile(bdata.ITI(vios|postvios),90)];
fprintf('ITI cutoffs (s): %.4f, %.4f\n', cutoff)
ITI = bdata.ITI;
ITI(ITI<min(cutoff) | ITI>max(cutoff)) = nan;

% time window for auc
T = linspace(-5, 10, 7229);
twin = [0 0.5]; % auc time window
[~,t0] = min(abs(T-twin(1)));
[~,t1] = min(abs(T-twin(2)));

usethese = vios;

diti_all = [-diff(ITI); nan];
diti = diti_all(usethese);

% bin diti into bins
nbin=6;
fprintf('slower: %i, faster: %i\n', sum(diti<0), sum(diti>0))
edges_neg = quantile(diti(diti<0), linspace(0, 1, nbin/2+1));
edges_pos = quantile(diti(diti>0), linspace(0, 1, nbin/2+1));
edges = [edges_neg(1:end-1), 0, edges_pos(2:end)];
disp(edges)
diti_bins = discretize(diti, edges);

DA = pdata(usethese,:);

auc_mean = nan(1,nbin);
figure; hold on
set(gcf, position=[473,269,419,338])
mycolors = [linspace(0.85, 0, nbin)', linspace(0.85, 0, nbin)', linspace(0.85, 0, nbin)'];
for bb=1:nbin
    DA_bc = baselineCorrect(T, -0.1, 0, DA(diti_bins==bb,:));
    plotPretty(T, DA_bc, mycolors(bb,:)); xlim([-0.1 0.8]); drawnow;
    auc = trapz(T(t0:t1), DA_bc(:,t0:t1)')';
    % filter out extreme outliers in auc
    thresh = [prctile(auc,0.1), prctile(auc,99.9)];
    auc(auc<min(thresh) | auc>max(thresh)) = nan;
    auc_mean(bb) = mean(auc, 'omitnan');
end

end


