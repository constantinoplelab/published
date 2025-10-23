function plot_fig3h_left(datadir)
% Plot z-scored dopamine at the reward cue conditioned on trial initiation
% time on the subsequent trial, averaged across rats. N = 10 rats.
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

nbin = 4;
fname = fullfile(datadir, 'data-published/DAbyNextITI_SOFF.mat');
if exist(fname, 'file')
    load(fname, 'da', 'ratList')
else
    % load all DA rats
    ratList = loadRats(datadir, 'da');
    
    da = cell(1,nbin);
    bins = nan(length(ratList),nbin);
    for rr=1:length(ratList)
        ratname = ratList{rr};
        fprintf('%i of %i: %s\n', rr, length(ratList), ratname)
        [da(rr,:), bins(rr,:)] = DAbyNextITI_SOFF(datadir, ratname, nbin);
    end
    save(fname, 'da', 'bins', 'ratList')
end


mycolors = {'#303030', '#555555', '#7a7a7a', '#9e9e9e'};
T = linspace(-5, 10, 7229);
twin = [-0.1 0.5]; 
[~,t0] = min(abs(T-twin(1)));
[~,t1] = min(abs(T-twin(2)));

figure
set(gcf, units='inches', renderer = 'painters', ...
    position=[9,4,2, 1.6])
for bb=1:nbin
    plotPretty(T(t0:t1), cell2mat(da(:,bb)), mycolors{bb});
    hold on
end
xlim([-0.1, 0.5])
xline(0, 'k--')
set(gca, tickdir='out', box='off', ytick=[0:0.2:0.4], fontsize=8)
xlabel('Time from Reward Cue (s)')
ylabel('DA')
subtitle(sprintf('N = %i rats', length(ratList)))

end

function [avgDA, avgX] = DAbyNextITI_SOFF(datadir, rat, nbins)

zscore = @(x) (x-mean(x,'omitnan'))/std(x,'omitnan');

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

% get [trial number, session, next ITI]
bsoff = bstruct.(rat).SideOff;
if iscell(bsoff.UniqueDay(1))
    bsoff.UniqueDay = str2double(bsoff.UniqueDay);
end
soff_trialinfo = [bsoff.TrialNumber, bsoff.UniqueDay];
% get next trial initiation time from CPIn behavior data and append to
% soff_trialinfo
bdata = bstruct.(rat).CPIn; 
if iscell(bdata.UniqueDay(1))
    bdata.UniqueDay = str2double(bdata.UniqueDay);
end
% nan ITIs on first trial
bdata.ITI(bdata.TrialNumber==1) = nan;

% determine iti cutoffs using all sessions
cutoff = [prctile(bdata.ITI,5), prctile(bdata.ITI,95)];
fprintf('ITI cutoffs (s): %.4f, %.4f\n', cutoff)

n = size(soff_trialinfo, 1);
match = nan(n, 1);  % Preallocate with NaN
for i = 1:n
    ix = find(bdata.TrialNumber == soff_trialinfo(i,1)+1 & ...
              bdata.UniqueDay == soff_trialinfo(i,2), 1, 'first');
    if ~isempty(ix)
        match(i) = ix;
    end
end
iti_values = nan(size(match));
iti_values(~isnan(match)) = bdata.ITI(match(~isnan(match)));
soff_trialinfo(:,3) = iti_values;

% time window for auc
T = linspace(-5, 10, 7229);
twin = [-0.1 0.5]; % auc time window
[~,t0] = min(abs(T-twin(1)));
[~,t1] = min(abs(T-twin(2)));

sessions = unique(bdata.UniqueDay);
xs_all = nan(length(sessions), nbins);
da_avg = cell(1,nbins);
for sess=1:length(sessions)
    % get bdata for this session
    these = soff_trialinfo(:,2) == sessions(sess);
    nextITI = soff_trialinfo(these,3);
    nextITI(nextITI<min(cutoff) | nextITI>max(cutoff)) = nan;
    nextITI = zscore(nextITI);

    % bin ITI
    xdata = nextITI;
    edges = quantile(xdata, [0:1/nbins:1]);
    bins = discretize(xdata, edges);    
    xs_all(sess,:) = (edges(1:end-1)+edges(2:end))/2;

    % get pdata for this session
    peh = pdata(bsoff.UniqueDay==sessions(sess),:);

    for bb=1:nbins
        da_bc = baselineCorrect(T, -0.1, 0, peh(bins==bb,:));
        da_avg{bb}(sess,:) = mean(da_bc(:,t0:t1), 1, 'omitnan');
    end
end

avgDA = cellfun(@(x) mean(x,1,'omitnan'), da_avg, UniformOutput=false);
avgX = mean(xs_all,1,'omitnan');

figure
set(gcf, units='inches', renderer='painters',...
    position=[6,4,2.65,2])
mycolors = {'#303030', '#555555', '#7a7a7a', '#9e9e9e'};
for bb=1:nbins
    plotPretty(T(t0:t1), da_avg{bb}, mycolors{bb});
    hold on
end
ylabel('Reward Cue DA')
xlabel('Time (s)')
xlim([-0.1 0.5])
set(gca, tickdir='out', box='off', fontsize=11)
xline(0, 'k--')
subtitle(rat)
drawnow

end


