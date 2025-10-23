function plot_fig3e_left(datadir)
% Plot z-scored dopamine at the offer cue conditioned on trial initiation
% time on the subsequent trial, averaged across rats. N = 10 rats.
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

nbin = 4;
fname = fullfile(datadir, 'data-published/DAbyNextITI.mat');
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
        [da(rr,:), bins(rr,:)] = DAbyNextITI(datadir, ratname, nbin);
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
set(gca, tickdir='out', box='off', ytick=[-0.2:0.2:0.2], fontsize=8)
xlabel('Time from Offer Cue (s)')
ylabel('DA')
subtitle(sprintf('N = %i rats', length(ratList)))

end

function [avgDA, avgX] = DAbyNextITI(datadir, rat, nbins)

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
if iscell(bdata.UniqueDay(1))
    bdata.UniqueDay = str2double(bdata.UniqueDay);
end
zscore = @(x) (x-mean(x,'omitnan'))/std(x,'omitnan');

% nan ITIs on first trial
bdata.ITI(bdata.TrialNumber==1) = nan;

% determine iti cutoffs using all sessions
cutoff = [prctile(bdata.ITI,5), prctile(bdata.ITI,95)];
fprintf('ITI cutoffs (s): %.4f, %.4f\n', cutoff)

T = linspace(-5, 10, 7229);
twin = [-0.1 0.5]; 
[~,t0] = min(abs(T-twin(1)));
[~,t1] = min(abs(T-twin(2)));

sessions = unique(bdata.UniqueDay);
xs_all = nan(length(sessions), nbins);
da_avg = cell(1,nbins);
for sess=1:length(sessions)
    % get bdata for this session
    beh = bdata(bdata.UniqueDay==sessions(sess),:);

    % DA on violation trial N vs ITI(N+1)
    usethese = beh.TrialType==2;
    
    % get next trial ITI
    iti = beh.ITI;
    iti(iti<min(cutoff) | iti>max(cutoff)) = nan;
    iti = zscore(iti);
    nextITI = [iti(2:end); nan];
    nextITI = nextITI(usethese);

    % bin ITI
    xdata = nextITI;
    edges = quantile(xdata, [0:1/nbins:1]);
    bins = discretize(xdata, edges);    
    xs_all(sess,:) = (edges(1:end-1)+edges(2:end))/2;

    % get pdata for this session
    peh = pdata(bdata.UniqueDay==sessions(sess),:);

    da = peh(usethese,:);
    for bb=1:nbins
        da_bc = baselineCorrect(T, -0.1, 0, da(bins==bb,:));
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
ylabel('Offer Cue DA')
xlabel('Time (s)')
xlim([-0.1 0.5])
set(gca, tickdir='out', box='off', fontsize=11)
xline(0, 'k--')
subtitle(rat)
drawnow

end


