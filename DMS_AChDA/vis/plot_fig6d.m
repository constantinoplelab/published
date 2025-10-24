function plot_fig6d(datadir)
% Plot side-on DA response on contralateral trials, bottom vs top quartile
% of reaction time distribution. N = 10 rats.
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

nbin = 4;
fname = fullfile(datadir, 'data-published/sideOnDAbyRT.mat');
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
        [da(rr,:), bins(rr,:)] = sideOnDAbyRT(datadir, ratname, nbin);
    end
    save(fname, 'da', 'ratList')
end

T = linspace(-5, 10, 7229);
twin = [0 0.5]; 
[~,t0] = min(abs(T-twin(1)));
[~,t1] = min(abs(T-twin(2)));
taxis = T(t0:t1);

M = cell(1,4);
for col=1:nbin
    M{col} = cell2mat(da(:,col));
end

% stats (signrank)
pvals = zeros(1, t1-t0+1);
for t=1:length(pvals)
    pvals(t) = signrank(M{1}(:,t), M{nbin}(:,t));
end
sig = pvals<0.05 & pvals>=0.01;

mycolors.fast = '#FF7252';
mycolors.slow = '#FFCAA6';

figure
set(gcf, units='inches', renderer = 'painters', ...
    position=[9,4,2, 1.6])

plotPretty(taxis, M{1}, mycolors.fast);
hold on
plotPretty(taxis, M{nbin}, mycolors.slow);

xlim(twin)
xline(0, 'k--')
yl = ylim;
ylim([yl(1) yl(2)+0.1])
plot(taxis(sig), yl(2)+0.05, 'k.')
set(gca, tickdir='out', box='off', fontsize=8)
xlabel('Reward Port On (s)')
ylabel('DA')
subtitle(sprintf('N = %i rats', length(ratList)))

end

function [avgDA, avgX] = sideOnDAbyRT(datadir, rat, nbin)

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
pdata = pstruct.(rat).SideOn;
load(fullfile(bfile.folder, bfile.name), 'bstruct')
bdata = bstruct.(rat).SideOn;

if iscell(bdata.UniqueDay(1))
    bdata.UniqueDay = str2double(bdata.UniqueDay);
end

T = linspace(-5, 10, 7229);
twin = [0 0.5]; 
[~,t0] = min(abs(T-twin(1)));
[~,t1] = min(abs(T-twin(2)));

sessions = unique(bdata.UniqueDay);
xs_all = nan(length(sessions), nbin);
da_avg = cell(1,nbin);

ch = extract_channel(pfile.name);

for sess=1:length(sessions)
    % get bdata, pdata for this session
    beh = bdata(bdata.UniqueDay==sessions(sess),:);
    peh = pdata(bdata.UniqueDay==sessions(sess),:);

    % get contra trials
    contra = beh.RewardPort~=ch;
    
    % display range of side led reaction time (s) if min or max is too
    % short/long
    if min(beh.SLRT(contra)) < 0.05 | max(beh.SLRT(contra)) > 20
        fprintf('SLRT (contra): %.3f ~ %.3f\n', ...
            min(beh.SLRT(contra)), max(beh.SLRT(contra)))
    end

    % bin SLRT using contralateral trials only
    xdata = beh.SLRT(contra);
    edges = quantile(xdata, [0:1/nbin:1]);
    bins = discretize(xdata, edges);    
    xs_all(sess,:) = (edges(1:end-1)+edges(2:end))/2;
    
    peh_contra = peh(contra,:);
    for bb=1:nbin
        da = peh_contra(bins==bb,:);
        da_bc = baselineCorrect(T, -0.1, 0, da);
        
        da_avg{bb}(sess,:) = mean(da_bc(:,t0:t1), 1, 'omitnan');
    end

end

avgDA = cellfun(@(x) mean(x,1,'omitnan'), da_avg, UniformOutput=false);
avgX = mean(xs_all,1,'omitnan');

end

function ch = extract_channel(filename)
    % Extract integer after 'ch' using regular expression
    tokens = regexp(filename, 'ch(\d+)', 'tokens');
    ch = str2double(tokens{1}{1});

end