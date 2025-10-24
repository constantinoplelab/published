function plot_fig5g(datadir)

fpath = fullfile(datadir, 'data-published/NpxData/pD1pD2.mat');
load(fpath, 'SU', 'S')

% example cells: 
% (session, neuron, trial) = (8,11,19), (10,2,16)
example = [8,11,19; 10,2,16]; 
mycolors = {'#B8B8B8', '#5B5B5B', 'k'};

figure;
t = tiledlayout(2,3);
set(gcf, units='inches', renderer='painters', position=[8,5,4,2])
for ex = 1:size(example,1)
    neuron = SU{example(ex,1)}{example(ex,2)};
    T = neuron.xvec.COFF;
    beh = S{example(ex,1)};
    [FR, validTrials, RPE_firstTrial] = ...
        get_FR_persistence(neuron, beh, 'COFF');
    data = FR{validTrials(example(ex,3))};
    for ii=1:3
        nexttile(ii+3*(ex-1))
        plot(T, data(ii,:), color=mycolors{ii}, linewidth=2);
        if ex==1
            ylim([0 20])
        else
            ylim([0 30])
        end
        if ii==1
            subtitle(sprintf('RPE=%.2f', RPE_firstTrial(example(ex,3))))
        end

        hold on
    end
    
end

for a=1:6
    nexttile(a);
    set(gca, box='  off', tickdir='out')
    xlim([-0.1 1])
    xline(0, 'k--')
end
xlabel(t, 'Offer Cue (s)', fontsize=9)
ylabel(t, 'FR (Hz)', fontsize=9)
drawnow;


end