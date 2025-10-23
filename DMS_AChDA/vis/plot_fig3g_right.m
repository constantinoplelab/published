function plot_fig3g_right(datadir)
% Plot z-scored rat trial initiation time and model simulation regressed
% against reward delays in the last 15 trials. For visualization, only 7
% trials are shown. N = 10 rats.
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

% regression settings
z_arg = 0;
nback=15; % number of trials back to use for regression
nback_vis = 7; % number of trials to show for visualization

%% Simulate trial initiation with a vanilla TD model
alpha = 0.46;
discount = 0.54;

[modelITI, ~, ratTrial, ~] = ...
    run_vanillaTD(datadir, 'delay', alpha, discount);
ratList = fields(ratTrial);

betas_model = nan(nback, length(ratList));

%% Regress model simulated trial initiation time against reward delay
for rr=1:length(ratList)
    ratname = ratList{rr};
    A = ratTrial.(ratname);
    A.ITI = modelITI.(ratname);

    [betas_model(:,rr),~] = regress_latency_vs_delay(A, nback, z_arg);
end

%% Regress rat ITI against reward delays
folder = fullfile(datadir, 'data-published/A_structs');
ratList = loadRats(datadir, 'da');
betas_rat = nan(nback, length(ratList));
for rr=1:length(ratList)
    rat = ratList{rr};
    file = strcat('ratTrial_', rat);
    load(fullfile(folder, file), 'A');
    [betas_rat(:,rr), ~] = regress_latency_vs_delay(A, nback, z_arg);
end

%% Plot
figure;
set(gcf, units='inches', renderer='painters', ...
    position=[6,4,2.8,3])

plotPretty(1:nback_vis, betas_rat(1:nback_vis,:)', '#793cab');
hold on
plotPretty(1:nback_vis, betas_model(1:nback_vis), 'k');
ylim([-0.25 0.1])
yticks([-0.2:0.1:0.1])
ylabel('Regression weight')
xlim([0.5 nback_vis+0.5])
xticks(1:nback_vis)
xlabel('Trials back')
yline(0, 'k--', 'linewidth', 0.75)
set(gca, box='off', tickdir='out', fontsize=11)
subtitle('Trial initiation time vs reward delay')


end

