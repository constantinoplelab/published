function plot_fig3g_left(datadir)
% Plot dopamine AUC at the reward cue (0-0.5 s) and model predicted RPE
% regressed against reward delays in the last 15 trials during mixed
% blocks. For visualization, only 7 trials back are shown. N = 10 rats.
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

% regression settings
nback=15; % number of trials back to use for regression
nback_vis = 7; % number of trials to show for visualization

%% Simulate RPE with a vanilla TD model
alpha = 0.46;
discount = 0.54;

[~, ~, ratTrial, modelRPE] = ...
    run_vanillaTD(datadir, 'delay', alpha, discount);

ratList = fields(modelRPE);

% Compute regression coefficients for model RPE
betas_model = nan(nback, length(ratList));
for rr=1:length(ratList)
    ratname = ratList{rr};
    A = ratTrial.(ratname);
    rpe_model = modelRPE.(ratname);
    % Regress model simulated RPE against reward history
    betas_model(:,rr) = regress_modelRPE_vs_delay(A, rpe_model, nback);
end

%% regress dopamine AUC against reward history
fpath = fullfile(datadir, 'data-published', 'DAAUC_regress_delay.mat');

if ~exist(fpath, 'file')
    % load list of DA rats
    ratList = loadRats(datadir, 'da');
    weights = cell(1,length(ratList));
    
    for r=1:length(ratList)
        rat = ratList{r};
        disp(rat)
        weights{r} = regressAUCvsDelay(datadir, rat, nback);
    end
    save(fpath, 'ratList', 'weights')
else
    load(fpath, 'ratList', 'weights');
end

% get session average of regression weights for each rat
coeff_avg = cellfun(@(x) mean(x,2,'omitnan'), weights, uniformoutput=false);
coeff_avg = cell2mat(coeff_avg);

%% Plot
figure;
set(gcf, units='inches', renderer='painters', ...
    position=[6,4,2.8,3])

plotPretty(0:nback_vis, coeff_avg(1:nback_vis+1,:)', '#793cab');
hold on
ylim([-0.01 0.03])
yticks([0, 0.02])
yline(0, 'k--', linewidth=0.75)
ylabel('Regression weight')

yyaxis right
plotPretty(0:nback_vis, betas_model(1:nback_vis+1,:)', 'k');
ylim([-0.35 1.05])
yticks([0:0.4:1])

ax = gca;
ax.YAxis(1).Color = '#793cab';
ax.YAxis(2).Color = 'k';
xlim([-0.5 nback_vis+0.5])
xticks(0:nback_vis)
xlabel('Trials back')
set(gca, box='off', tickdir='out', fontsize=11)
subtitle('RPE vs reward delay')

end

function betas = regressAUCvsDelay(datadir, rat, nback, mixedblock)

if nargin<4
    mixedblock = 0; % use all blocks by default
end

T = linspace(-5, 10, 7229);
[~,t0] = min(abs(T));
[~,t1] = min(abs(T-0.5));

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
pdata = pstruct.(rat).SideOff;
load(fullfile(bfile.folder, bfile.name), 'bstruct')
bdata = bstruct.(rat).SideOff;

% regress each session
newsession = [1; find(diff(bdata.TrialNumber)<0)+1; size(pdata,1)+1];
betas = nan(nback+2, length(newsession)-1);

for ss=1:length(newsession)-1
    sessionP = pdata(newsession(ss):newsession(ss+1)-1, :);
    sessionB = bdata(newsession(ss):newsession(ss+1)-1,:);
    % baseline correct
    sessionP = baselineCorrect(T, -0.1, 0, sessionP);
    sessionP = sessionP(:,t0:t1);
    auc = trapz(T(t0:t1), sessionP, 2);
            
    % get reward history matrix
    RD = sessionB.RewardDelay;
    RD(RD==100) = nan;
    
    if mixedblock==1
        auc(sessionB.Block~=1) = nan;
        RD(sessionB.Block~=1) = nan;
    end

    % Create design matrix
    X = nan(length(RD), nback+2);
    for j = 1:nback+1
        X(:,j) = [zeros(j-1,1); RD(1:end-j+1)];    
    end
    X(:,end) = ones(length(RD),1); %constant term.

    betas(:,ss) = regress(auc, X);

end

end