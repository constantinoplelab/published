function plot_fig3d_left(datadir)
% Plot dopamine AUC at the offer cue (0-0.5 s) and model predicted RPE
% regressed against reward offers in the last 15 trials during mixed
% blocks. For visualization, only 7 trials back are shown. N = 10 rats.
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

% regression settings
hit_arg = 0; % don't zero-out unrewarded trials
allblocks = 0; % mixed block only
nback=15; % number of trials back to use for regression
nback_vis = 7; % number of trials to show for visualization

%% Simulate RPE with a vanilla TD model
alpha = 0.5;
discount = 0.48;

[~, ~, ratTrial, modelRPE] = ...
    run_vanillaTD(datadir, 'reward', alpha, discount);

ratList = fields(modelRPE);

% Compute regression coefficients for model RPE
betas_model = nan(nback+2, length(ratList));
for rr=1:length(ratList)
    ratname = ratList{rr};
    A = ratTrial.(ratname);
    rpe_model = modelRPE.(ratname);
    % Regress model simulated RPE against reward history
    betas_model(:,rr) = regress_modelRPE_vs_rew(A, rpe_model, nback, ...
        hit_arg, allblocks);
end

%% regress dopamine AUC against reward history
fpath = fullfile(datadir, 'data-published', 'DAAUC_regress_rew.mat');

if ~exist(fpath, 'file')
    % load list of DA rats
    ratList = loadRats(datadir, 'da');
    weights = cell(1,length(ratList));
    
    for r=1:length(ratList)
        rat = ratList{r};
        disp(rat)
        weights{r} = regressAUCvsReward(datadir, rat, nback, hit_arg, allblocks);
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
ylim([-0.025 0.075])
yticks([-0.03:0.03:0.08])
ylabel('Regression weight')

yyaxis right
plotPretty(0:nback_vis, betas_model(1:nback_vis+1,:)', 'k');
ylim([-1.2 3.5])
yticks([-1:1:4])

ax = gca;
ax.YAxis(1).Color = '#793cab';
ax.YAxis(2).Color = 'k';
xlim([-0.5 nback_vis+0.5])
xticks(0:nback_vis)
xlabel('Trials back')
yline(0, 'k--', 'linewidth', 0.75)
set(gca, box='off', tickdir='out', fontsize=11)
subtitle('RPE vs reward history')

end

function betas = regressAUCvsReward(datadir, rat, nback, hit_arg, allblocks)

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
pdata = pstruct.(rat).CPIn;
load(fullfile(bfile.folder, bfile.name), 'bstruct')
bdata = bstruct.(rat).CPIn;

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
    
    if allblocks==0
        auc(sessionB.Block~=1) = nan;
    end
        
    % get reward history matrix
    sessionB.Reward = log2(sessionB.Reward);
    if hit_arg==1 % zero-out non-hit rewards
        hit = sessionB.TrialType==1;
        R = sessionB.Reward.*hit;
    else
        R = sessionB.Reward;
    end
    
    if allblocks==0
        R(sessionB.Block~=1) = nan;
    end

    % Create design matrix
    X = nan(length(R), nback+2);
    for j = 1:nback+1
        X(:,j) = [zeros(j-1,1); R(1:end-j+1)];    
    end
    X(:,end) = ones(length(R),1); %constant term.

    betas(:,ss) = regress(auc, X);
end

end