function plot_fig6f(datadir)
% Plot opt-out DA response on contralateral trials, bottom vs top quartile
% of next trial ITI distribution. N = 10 rats.
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

nbin = 4;
fname = fullfile(datadir, 'data-published/optOutDAbyNextITI.mat');
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
        [da(rr,:), bins(rr,:)] = optOutDAbyNextITI(datadir, ratname, nbin);
    end
    save(fname, 'da', 'ratList')
end


T = linspace(-5, 10, 7229);
twin = [-0.3 0.2]; 
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
sig = pvals<0.05;


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
ylim([yl(1) 0.4])
yticks(0:0.2:0.4)
plot(taxis(sig), 0.38, 'k.')
set(gca, tickdir='out', box='off', fontsize=8)
xlabel('Opt Out (s)')
ylabel('DA')
subtitle(sprintf('N = %i rats', length(ratList)))


end

function [avgDA, avgX] = optOutDAbyNextITI(datadir, rat, nbin)

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
pdata = pstruct.(rat).OptOut;
load(fullfile(bfile.folder, bfile.name), 'bstruct')

% get [trial number, session, next ITI]
beh_opt = bstruct.(rat).OptOut;
if iscell(beh_opt.UniqueDay(1))
    beh_opt.UniqueDay = str2double(beh_opt.UniqueDay);
end
opt_trialinfo = [beh_opt.TrialNumber, beh_opt.UniqueDay];
% get next trial initiation time from CPIn behavior data and append to
% opt_trialinfo
beh_coff = bstruct.(rat).CPIn; 
if iscell(beh_coff.UniqueDay(1))
    beh_coff.UniqueDay = str2double(beh_coff.UniqueDay);
end
% nan ITIs on first trial
beh_coff.ITI(beh_coff.TrialNumber==1) = nan;

% determine iti cutoffs using all sessions
cutoff = [prctile(beh_coff.ITI,1), prctile(beh_coff.ITI,98)];
fprintf('ITI cutoffs (s): %.4f, %.4f\n', cutoff)

n = size(opt_trialinfo, 1);
match = nan(n, 1);  % Preallocate with NaN
for i = 1:n
    ix = find(beh_coff.TrialNumber == opt_trialinfo(i,1)+1 & ...
              beh_coff.UniqueDay == opt_trialinfo(i,2), 1, 'first');
    if ~isempty(ix)
        match(i) = ix;
    end
end
iti_values = nan(size(match));
iti_values(~isnan(match)) = beh_coff.ITI(match(~isnan(match)));
iti_values(iti_values<min(cutoff) | iti_values>max(cutoff)) = nan;
opt_trialinfo(:,3) = iti_values;

ch = extract_channel(pfile.name);

% use session aggregate data to bin trials since opt-out trials are limited
contra = beh_opt.RewardPort~=ch;
xdata = opt_trialinfo(contra,3);
edges = quantile(xdata, [0:1/nbin:1]);
bins = discretize(xdata, edges);    
avgX = (edges(1:end-1)+edges(2:end))/2;

T = linspace(-5, 10, 7229);
twin = [-0.3 0.2]; 
[~,t0] = min(abs(T-twin(1)));
[~,t1] = min(abs(T-twin(2)));

peh_contra = pdata(contra, :);
for bb=1:nbin
    da = peh_contra(bins==bb,:);
    da_bc = baselineCorrect(T, -0.3, 0, da);

    avgDA{1,bb} = mean(da_bc(:,t0:t1), 1, 'omitnan');   
end

figure
set(gcf, units='inches', renderer='painters',...
    position=[6,4,2.65,2])
mycolors.fast = '#FF7252';
mycolors.slow = '#FFCAA6';

plot(T(t0:t1), avgDA{1}, color=mycolors.fast);
hold on
plot(T(t0:t1), avgDA{nbin}, color=mycolors.slow);

ylabel('DA')
xlabel('Opt Out (s)')
xlim(twin)
set(gca, tickdir='out', box='off', fontsize=11)
xline(0, 'k--')
subtitle(rat)
drawnow

end

function ch = extract_channel(filename)
    % Extract integer after 'ch' using regular expression
    tokens = regexp(filename, 'ch(\d+)', 'tokens');
    ch = str2double(tokens{1}{1});

end