function plot_fig5h(datadir)

prct = 15;
fpath = fullfile(datadir, 'data-published/NpxData/pD1pD2.mat');
load(fpath, 'SU', 'S', 'pd1_coff_list', 'pd2_coff_list');
plot_persistence(SU, S, pd1_coff_list, pd2_coff_list, prct)

end

function plot_persistence(SU, S, d1list, d2list, prct)

if nargin<5
    prct = 15;
end
AvgFRs = struct;
RPEinit = struct;
AvgFRs_norm = struct;
twin = [0 0.3];

% D1
[AvgFRs.D1, AvgFRs_norm.D1, RPEinit.D1, ~] = ...
    getAllData(SU, S, d1list, twin);

% D2
[AvgFRs.D2, AvgFRs_norm.D2, RPEinit.D2, ~] = ...
    getAllData(SU, S, d2list, twin);

%%
figure;
set(gcf, units='inches', render='painters', ...
    position=[8,5,2,2.7])
msn = fields(RPEinit);
mycolors = {'#004489', '#b92237'};

% Positive RPE on trial N-1
p = nan(2,3);
n = [];
for m=1:length(msn)
    ub = prctile(RPEinit.(msn{m}), 100-prct);
    these = find(RPEinit.(msn{m})>=ub & ~isnan(AvgFRs_norm.(msn{m})(:,1)));
    n(m) = length(these);
    errorbar(1:3, mean(AvgFRs_norm.(msn{m})(these,:),1,'omitnan'), ...
        sem(AvgFRs_norm.(msn{m})(these,:),'omitnan'), '-', ...
        color=mycolors{m}, capsize=0, linewidth=1.5);
    hold on
    p(m,1) = signrank(AvgFRs.(msn{m})(these,1), AvgFRs.(msn{m})(these,2));
    p(m,2) = signrank(AvgFRs.(msn{m})(these,1), AvgFRs.(msn{m})(these,3));
    p(m,3) = signrank(AvgFRs.(msn{m})(these,2), AvgFRs.(msn{m})(these,3));

end
disp(p)

set(gca, box='off', tickdir='out', fontsize=13, xtick=[1,2,3], ...
    xticklabels={'N-1', 'N', 'N+1'})
xlim([0.75 3.25])
xlabel('Trial index')
ylabel('Average FR (norm.)')

    
end

function [AvgFRs, AvgFRs_norm, RPEinit, timeDur] = ...
    getAllData(SU, S, sigCells, twin)

AvgFRs = cell(length(sigCells),1);
RPEinit = cell(length(sigCells),1); % RPE on trial N-1
timeDur = cell(length(sigCells),1); % time from COFF on trial N-1 to N+1 (s)
for cc=1:length(sigCells)
    example = sigCells(cc,:);
    neuron = SU{example(1)}{example(2)};
    beh = S{example(1)};
    [avgFR, RPE_firstTrial, trialIdx] = ...
        get_avgFR_persistence(neuron, beh, 'COFF', twin);
    AvgFRs{cc} = avgFR;
    RPEinit{cc} = RPE_firstTrial;
    timeDur{cc} = arrayfun(@(xx) beh.Cled(xx+2,2)-beh.Cled(xx,2), trialIdx);
end
RPEinit = cell2mat(RPEinit);
AvgFRs = cell2mat(AvgFRs);
timeDur = cell2mat(timeDur);
AvgFRs_norm = AvgFRs ./ AvgFRs(:,1);


end