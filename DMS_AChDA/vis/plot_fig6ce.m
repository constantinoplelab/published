function plot_fig6ce(datadir)
% Plot difference in time to peak in dopamine compared to head speed (N =
% 24 sessions from 5 rats) when the side LED turns on, and when rats move
% towards the opt-out port to opt-out.
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.


% Load rats with DLC recordings
rats = {'M036', 'J027', 'J029', 'J063', 'J059'};
load(fullfile(datadir, 'data-published/ratlist.mat'), 'channelList');

A = {'SideOn', 'OptOut'};
latV = cell(length(rats),1);
latF = cell(length(rats),1);

for r=1:length(rats)
    ratname = rats{r};
    channel = channelList.(ratname);
    [latV{r}, latF{r}] = ...
        getLatency_speed_photometry(datadir, ratname, channel, A, false);
end

% Plot 
latF = cell2mat(latF);
latV = cell2mat(latV);

N = size(latF,1);
jitx = rand(N,1)*0.1;
jitx(1:N/2) = -jitx(1:N/2);

figure;
set(gcf, units='inches', renderer='painters', position=[8,4,3,2])
t = tiledlayout(1,length(A));
p = nan(1,length(A));
for a=1:length(A)
    nexttile(a); hold on
    plot(1+jitx, latV(:,a)-latF(:,a), 'o', color='#a0a0a0');
    errorbar(1, mean(latV(:,a)-latF(:,a)), sem(latV(:,a)-latF(:,a)), ...
        'k_', capsize=0, linewidth=1.5)
    yline(0, 'k--')
    xlim([0.75 1.25])
    ylim([-410 410])
    set(gca, xtick=[], fontsize=8, tickdir='out', box='off')
    [~,p(a)] = ttest(latV(:,a)-latF(:,a));
    subtitle(sprintf('p=%.3f', p(a)))
end
ylabel(t, 'Speed - DA (ms)', fontsize=8)
nexttile(1); xlabel('Reward Port On')
nexttile(2); xlabel('Opt Out')

end

function [latV, latF] = ...
    getLatency_speed_photometry(datadir, ratName, channel, A, plotarg)

if nargin<5
    plotarg=false;
end

files = dir(fullfile(datadir, 'data-published/DLCdata', ...
    strcat(ratName, '_*_DLC.mat')));

latV = nan(length(files), length(A));
latF = nan(length(files), length(A));

for f=1:length(files)
    % load data
    date = regexp(files(f).name, '\d{8}', 'match');
    date = date{1};
    [DLC, pData, bData] = loadAllStructs(datadir, ratName, date);
    pData = pData.AlignmentStruct.Photometry;

    tWin_dlc = [0 0.7];
    T_dlc = DLC.Times;
    [~,t0dlc] = min(abs(T_dlc-tWin_dlc(1)));
    [~,t1dlc] = min(abs(T_dlc-tWin_dlc(2)));

    T_p = pData.Times;

    if plotarg
        figure;
        t = tiledlayout(1,length(A));
        set(gcf, render='painters', units='inches', ...
            position=[7,5.5,5.48,2.97])
    end
    for a=1:length(A)
        if strcmpi(A{a}, 'OptOut')
            tWin_p = [-0.3 0.4];
        else
            tWin_p = [0 0.7];
        end
        [~,t0p] = min(abs(T_p-tWin_p(1)));
        [~,t1p] = min(abs(T_p-tWin_p(2)));

        % get contra trials
        if channel==1
            ix = bData.(A{a}).RewardPort ==2;
        else
            ix = bData.(A{a}).RewardPort ==1;
        end
        v = DLC.speed.(A{a})(ix, t0dlc:t1dlc);
        fl = pData.(A{a})(ix, t0p:t1p);

        % average
        v = mean(v, 1, 'omitnan');
        fl = mean(fl, 1, 'omitnan');
        if plotarg
            nexttile(a);
            plot(T_p, mean(pData.(A{a})(ix,:),1,'omitnan'), ...
                color='#ffa052', linewidth=2);
            
            xline(0, 'k--')
            hold on
            yyaxis right
            plot(T_dlc, mean(DLC.speed.(A{a})(ix,:),1,'omitnan'), ...
                color='k', linewidth=2);
            ax = gca;
            ax.YAxis(1).Color = '#ffa052';
            ax.YAxis(2).Color = 'k';
            subtitle(A{a})
            set(gca, box='off', tickdir='out')
            if strcmpi(A{a}, 'OptOut')
                xlim([-0.3 0.7]);
            else
                xlim([-0.1 0.9]);
            end
        end

        % get latency to peak
        [~,latV(f,a)] = max(v);
        latV(f,a) = T_dlc(latV(f,a)+t0dlc)*1000; 
        [~,latF(f,a)] = max(fl);
        latF(f,a) = T_p(latF(f,a)+t0p)*1000;
    end
    if plotarg
        title(t, sprintf('%s: %s', ratName, date))
        xlabel(t, 'Time (s)')
        drawnow;
    end
end

end

function [DLCdata, pdata, bData] = ...
    loadAllStructs(datadir, ratName, date)
    
    fprintf('%s: %s\n', ratName, date)

    if ~contains(date, '/')
        date = strcat(date(1:4),'/',date(5:6),'/',date(7:end));
    end

    folder = fullfile(datadir, 'data-published/DLCdata');

    % DLC file
    searchstr = strcat(ratName, '_*', datestr(date, 'yyyymmdd'), '_*_DLC.mat');
    file = dir(fullfile(folder, searchstr));
    file = file.name;
    DLCdata = load(fullfile(folder, file), 'AlignmentStruct');
    DLCdata = DLCdata.AlignmentStruct;
    disp('Loaded DLC data')

    % bData file
    searchstr = strcat(ratName, '_*', datestr(date, 'yyyymmdd'), '_*_bData.mat');
    file = dir(fullfile(folder, searchstr));
    file = file.name;
    load(fullfile(folder, file), 'bData');
    disp('Loaded bData')

    % photometry file
    searchstr = strcat(ratName, '_*', datestr(date, 'yyyymmdd'), '*_tmac.mat');
    file = dir(fullfile(folder, searchstr));
    file = file.name;
    pdata = load(fullfile(folder, file), 'AlignmentStruct');
    disp('Loaded photometry data')

    
end

