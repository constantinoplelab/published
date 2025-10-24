function plot_fig5c(datadir)
% Plot firing rate of two example cells on current and subsequent trial
% following a positive and negative RPE. 
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

% load data
fpath = fullfile(datadir, 'data-published/NpxData/pD1pD2.mat');
load(fpath, 'goodD1', 'SU', 'S')
neuronlist = goodD1;
dFRposRPE = cell(length(neuronlist),1);
dFRnegRPE = dFRposRPE;
dFRposRPE_SON = dFRposRPE;
dFRnegRPE_SON = dFRposRPE;
zscorefr = false;

for neuron=1:length(neuronlist)
    sess = neuronlist(neuron,1);
    s = neuronlist(neuron,2);
    [dFRpos, dFRneg, dFRposSON, dFRnegSON, tt] = ...
        get_dFRbyRPE(SU{sess}{s}, S{sess}, zscorefr);
    dFRposRPE{neuron} = dFRpos;
    dFRnegRPE{neuron} = dFRneg;
    dFRposRPE_SON{neuron} = dFRposSON;
    dFRnegRPE_SON{neuron} = dFRnegSON;
end

% example cells: (neuron, trial) 
dFRposRPE_good = [23,29; 22,46];
dFRnegRPE_good = [23,42; 22,52];

% plot: current vs next trial on positive/negative RPE trials for each neuron
for ex=1:length(dFRposRPE_good)
    neuron = dFRposRPE_good(ex,1);
    trial = dFRposRPE_good(ex,2);
    fig = plot_TrialsetFR(dFRposRPE{neuron}{trial}, ...
        dFRposRPE_SON{neuron}{trial}, tt);
    title(fig, {sprintf('Example DMS cell %i', ex), 'Positive RPE trial'}, ...
        fontsize=9)
    drawnow;
end
%
for ex=1:length(dFRnegRPE_good)
    neuron = dFRnegRPE_good(ex,1);
    trial = dFRnegRPE_good(ex,2);
    fig = plot_TrialsetFR(dFRnegRPE{neuron}{trial}, ...
        dFRnegRPE_SON{neuron}{trial}, tt);
    title(fig, {sprintf('Example DMS cell %i', ex), 'Negative RPE trial'}, ...
        fontsize=9)
    drawnow;
end

end

function [dFRposRPE, dFRnegRPE, dFRposRPE_SON, dFRnegRPE_SON, tt] = ...
    get_dFRbyRPE(neuron, S, zscorefr)
% get delta FR for big positive/negative RPE trials

if nargin<3
    zscorefr = false;
end

if zscorefr
    % z-score hmat
    rawFR = neuron.hmat.COFF;
    neuron.hmat.COFF = zscore(rawFR(:), 0, 'omitnan');
    neuron.hmat.COFF = reshape(neuron.hmat.COFF, size(rawFR));
end

A = {'COFF', 'SON'};

% plot twin sec around the event
twin = [-0.3 1.0];
T = neuron.xvec.COFF;
[~,t1] = min(abs(T-twin(1)));
[~,t2] = min(abs(T-twin(2)));

[~, d_iti, idx] = get_dFRdITI(neuron, S, A);

spikes = neuron.hmat.COFF(:, t1:t2);
spikes_SON = neuron.hmat.SON(:, t1:t2);

% find big positive RPE trials (last N bins)
N = 3;

xdata = d_iti{1};
nbins = 10;
edges = linspace(min(xdata), max(xdata), nbins+1);
bins = discretize(d_iti{1}, edges);

posRPE = find(bins>=nbins-N+1);

dFRposRPE = cell(length(posRPE),1);
dFRposRPE_SON = cell(length(posRPE),1);
for ii=1:length(posRPE)
    % e.g., posRPE(ii)=10, meaning d_iti(10)= -iti(11)+iti(10) is big.
    % Trial indices for these trials are m=idx(10) and n=idx(11).
    % For COFF, compare spikes on trial m vs n 
    % For SON, compare spikes on trial m vs the next time SON occurs (i.e.,
    % next non-violation trial)
    
    % COFF:
    m = idx{1}(posRPE(ii));
    n = idx{1}(posRPE(ii)+1);
    dFRposRPE{ii} = spikes([m,n],:);
    % SON:
    if ~isempty(find(idx{2}==m))
        if find(idx{2}==m) == length(idx{2})
            dFRposRPE_SON{ii}(1,:) = spikes_SON(m,:);
            dFRposRPE_SON{ii}(2,:) = nan;
        else
            n = idx{2}(find(idx{2}==m)+1);
            dFRposRPE_SON{ii} = spikes_SON([m,n],:);
        end
    else
        dFRposRPE_SON{ii} = nan(2,t2-t1+1);
    end
end

negRPE = find(bins<=N);
dFRnegRPE = cell(length(negRPE),1);
dFRnegRPE_SON = cell(length(negRPE),1);
for ii=1:length(negRPE)
    % COFF:
    m = idx{1}(negRPE(ii));
    n = idx{1}(negRPE(ii)+1);
    dFRnegRPE{ii} = spikes([m,n],:);
    % SON:
    if ~isempty(find(idx{2}==m))
        if find(idx{2}==m) == length(idx{2})
            dFRnegRPE_SON{ii}(1,:) = spikes_SON(m,:);
            dFRnegRPE_SON{ii}(2,:) = nan;
        else
            n = idx{2}(find(idx{2}==m)+1);
            dFRnegRPE_SON{ii} = spikes_SON([m,n],:);
        end
    else
        dFRnegRPE_SON{ii} = nan(2,t2-t1+1);
    end
end

tt = T(t1:t2);

end

function t = plot_TrialsetFR(dFRCOFF, dFRSON, tt)
% plot current and next FR at COFF and SON

mycolor = {'#999999', 'k'};


figure;
t = tiledlayout(1,2, Padding='compact');
set(gcf, units='inches', renderer='painters', position=[7,5,3.6,1.6])
yrange = [0 0];
nexttile(1);
for ii=1:2
    plot(tt, dFRCOFF(ii,:), color=mycolor{ii});
    hold on
end
xlabel('Offer Cue (s)')
yl = ylim;
yrange = updateYlim(yrange, yl);
nexttile(2);
for ii=1:2
    plot(tt, dFRSON(ii,:), color=mycolor{ii});
    hold on
end
xlabel('Reward Port On (s)')
yl = ylim;
yrange = updateYlim(yrange, yl);

for a=1:2
    nexttile(a);
    ylim(yrange);
    xlim([min(tt) max(tt)])
    xline(0, 'k--')
    xline(0.3, 'k--')
    set(gca, tickdir='out', box='off', fontsize=8)
end

ylabel(t, 'FR (Hz)', fontsize=8)

end











