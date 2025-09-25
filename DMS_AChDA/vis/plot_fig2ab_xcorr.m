function plot_fig2ab_xcorr(datadir)
% Plot cross-correlation of rat-averaged dopamine and acetylcholine
% signals at the offer cue on 40, 80uL trials in mixed block (left) and at
% the reward cue on different reward delay quartiles (right). 
% N = 10 DA rats, N = 10 ACh rats
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

A = {'Offer Cue', 'Reward Cue'}; 
region = 'dms';
datatypes = {'volume', 'delay'};

T = linspace(-5, 10, 7229);
dT = (T(3)-T(2))*1000; % ms
[~,t1] = min(abs(T+0.1));
[~,t2] = min(abs(T-0.5));
maxLag = round(500/dT);

figure;
set(gcf, renderer='painters', units='inches', ...
    position=[6,4,4,1.76]);
t = tiledlayout(1,2);

for d=1:2
    datatype = datatypes{d};
    [da_avg, ~] = get_avgFl(datadir, region, 'da', datatype);
    [ach_avg, ~] = get_avgFl(datadir, region, 'ach', datatype);
    
    mycolor = getColorScheme(datatype);
    if strcmp(datatype, 'volume')
        % get 40,80uL trials (4,5th rows) at offer cue event (2nd column)
        da_avg = da_avg(4:5,2);
        ach_avg = ach_avg(4:5,2);
        mycolor = mycolor(4:5);
    elseif strcmp(datatype, 'delay')
        % get reward cue event (4th column)
        da_avg = da_avg(:,4);
        ach_avg = ach_avg(:,4);
    end

    da_avg_trunc = cellfun(@(x) x(t1:t2), da_avg, UniformOutput=false);
    ach_avg_trunc = cellfun(@(x) x(t1:t2), ach_avg, UniformOutput=false);
    
    nexttile(d); hold on
    for r=1:size(da_avg_trunc,1)
        x = da_avg_trunc{r} - da_avg_trunc{r}(1);
        y = ach_avg_trunc{r} - ach_avg_trunc{r}(1);
       
        [cc, ll] = xcorr(y, x, maxLag, 'normalized');
        ll = ll*dT;
        plot(ll, cc, color=mycolor{r});
    end
    ylim([-1 1])
    xlim([-250 250])
    yline(0, 'k--')
    xlabel(A{d})
    xline(0, 'k-')
    set(gca, tickdir='out', box='off')
    
end

xlabel(t, 'Lag (ms)', fontsize=9)
ylabel(t, 'DA-ACh cross-corr.', fontsize=9)

end
