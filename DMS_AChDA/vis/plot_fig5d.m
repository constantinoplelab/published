function plot_fig5d(datadir)
% Plot trial-by-trial changes in average firing rate surrounding each task
% event (0-0.3 s) vs -delta trial initiaton time for an example cell.
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

fpath = fullfile(datadir, 'data-published/NpxData/pD1pD2.mat');
load(fpath, 'S', 'SU', 'offset', 'slope', 'pval')
session = 6;
neuron = 1;

n = cellfun(@(x) length(x), SU);
n = [0, cumsum(n)];

nIdx = n(session) + neuron;
coeff = [offset(:,nIdx), slope(:,nIdx)];
plot_dFR_dITI(SU{session}{neuron}, S{session}, coeff, pval(:,nIdx), 10)

end

function plot_dFR_dITI(neuron, S, coeff, pvals, nbins)

if nargin<5
    nbins = 7;
    if nargin<4
        pvals = [nan,nan];
    end
end

A = {'COFF', 'SON', 'SOFF'};
[d_fr, d_iti, idx] = get_dFRdITI(neuron, S, A);

xdata = d_iti{1};
edges = linspace(min(xdata), max(xdata), nbins+1);
midX = (edges(1:end-1)+edges(2:end))./2;
avgdFR = nan(nbins, length(A));
err = avgdFR;
for a=1:length(A)
    bin = discretize(d_iti{a}, edges);
    avgdFR(:,a) = arrayfun(@(bb) mean(d_fr{a}(bin==bb), 'omitnan'), ...
        1:nbins);
    err(:,a) = arrayfun(@(bb) sem(d_fr{a}(bin==bb), 'omitnan'), ...
        1:nbins);
end

% scatter plot
figure;
t = tiledlayout(2,length(A), padding='compact');
set(gcf, units='inches', position=[7,3.75,5.3,3.9], ...
    render='painters')
range = [0 0];
for a=1:length(A)
    nexttile;
    scatter(midX, avgdFR(:,a), MarkerEdgeColor='k'); 
    hold on
    errorbar(midX, avgdFR(:,a), err(:,a), capsize=0, ...
        color='k', linestyle='none');
    plot(midX, coeff(a,2)*midX+coeff(a,1), 'k');
    subtitle({A{a}, sprintf('%.3f (p=%.3f)', coeff(a,2), pvals(a))})
    set(gca, tickdir='out', box='off', fontsize=11)
    xlim([edges(1)-0.5 edges(end)+0.5])
    yl = ylim;
    range = updateYlim(range, yl);
end
for a=1:length(A)
    nexttile(a);
    ylim(range)
end
nexttile(2); 
xlabel('-\Delta Trial initiation time (s, log_{2})')

% big positive vs negative RPE trials
twin = [0 0.3];
T = neuron.xvec.(A{1});
[~,t1] = min(abs(T-twin(1)));
[~,t2] = min(abs(T-twin(2)));

N = 2;
range = [0 0];
for a=1:length(A)
    bin = discretize(d_iti{a}, edges);
    spikes = neuron.hmat.(A{a})(idx{a}, t1:t2);
    these = find(bin<=N);
    nexttile;
    plotPretty(T(t1:t2), spikes(these+1,:)-spikes(these,:), '#a0a0a0');
    hold on
    these = find(bin>=nbins-N+1);
    plotPretty(T(t1:t2), spikes(these+1,:)-spikes(these,:), 'k');
    yl = ylim;
    range = updateYlim(range, yl);
    set(gca, fontsize=11)
end
for a=length(A)+1:length(A)*2
    nexttile(a)
    ylim(range);
    yline(0, 'k--')
end
nexttile(5); 
xlabel('Time from event (s)')

ylabel(t, '\Delta FR (Hz)')
drawnow

end

