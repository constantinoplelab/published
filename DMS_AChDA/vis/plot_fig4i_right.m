function plot_fig4i_right(datadir)
% Plot z-scored latency to poke in the reward port after the reward port
% light turns on following short vs long reard delay trials, averaged
% across rats. N = 16 rats.
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

load(fullfile(datadir, 'data-published/ratlist.mat'), 'ratList');
folder = fullfile(datadir, 'data-published/A_Structs');

conditions = {'post short', 'post long'};
pokeLatency = nan(length(ratList), length(conditions)); 
for rr=1:length(ratList)
    rat = ratList{rr};
    fprintf('%i of %i: %s\n', rr, length(ratList), rat)
    file = strcat('ratTrial_', rat);
    load(fullfile(folder, file), 'A');

    bins = get_delayQuartile(datadir, rat);
    lb = bins(2); ub = bins(end-1);

    rd_prev = [nan; A.reward_delay(1:end-1)];
    rd_prev(rd_prev==100) = nan;
    posthit = [0; A.hits(1:end-1)];
    [A.reward,~] = convertreward(A.reward);

    mixed_posthit_nonvios = A.block==1 & posthit & ~A.vios;
    
    % filter outliers in side light reaction time & z-score
    cutoff = prctile(A.slrt, 99.99);
    fprintf('Filter out slrt > %.4f s\n', cutoff)
    A.slrt(A.slrt>cutoff) = nan;
    % z-score side light reaction time only using mixed block trials
    A.slrt(A.block~=1) = nan;
    A.slrt = (A.slrt - mean(A.slrt, 'omitnan')) ./ std(A.slrt, 'omitnan');
    for c=1:length(conditions)
        if contains(conditions{c}, 'short')
            these = mixed_posthit_nonvios & rd_prev <=lb;
        elseif contains(conditions{c}, 'long')
            these = mixed_posthit_nonvios & rd_prev >=ub;
        end
        pokeLatency(rr,c) = mean(A.slrt(these), 'omitnan');
    end
end

figure;
set(gcf, units='inches', position=[7,3.5,1.7,1.7])
plotPairs(pokeLatency(:,1), pokeLatency(:,2))
set(gca, xtick=[1 2], xticklabels=conditions, fontsize=8)
ylabel('Poke latency (z-score)')
xlim([0.5 2.5])
ylim([min(pokeLatency(:))-0.02 max(pokeLatency(:))+0.02])
pval = signrank(pokeLatency(:,1), pokeLatency(:,2));
disp(pval)

end
