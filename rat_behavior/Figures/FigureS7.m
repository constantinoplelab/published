function FigureS7(datadir, codedir)
%FigureS6 -  Alternative retrospective models fail to capture both fast 
% and slow trial initiation time dynamics at block transitions.
%   datadir = directory of dataset
%   codedir = director of code

%% Set up paths
s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codedir, 'IgnoreCase', ispc);

if ~onPath % only add code dir to path if it already isn't
    addpath(genpath(codedir))
end

%% Simulate data

% Simulation parameters
Nsims = 50; % number of simulations
twin = 40; % number of trials around each block to pull
smoothfactor = 20; % size of smoothing window
alphas = [0.05, 0.4, 0.8]; % range of learning rates

fname = 'ratTrial_J004.mat';
A = load([datadir 'A_Structs_Final' filesep fname]);
A = A.A;

[ltom_vanilla, htom_vanilla, ltom_rpe, htom_rpe] =...
    deal(cell(1, length(alphas)));

for aa = 1:length(alphas)
    alpha0 = alphas(aa);

    [A_mf_vanilla, A_mf_rpe] = deal(A);
    params = [alpha0, 0.2527];

    [ltom_vanilla{aa}, htom_vanilla{aa},...
        ltom_rpe{aa}, htom_rpe{aa}] = deal(nan(Nsims, 2*twin+1));

    rng(724)
    for rr = 1:Nsims
        disp([aa rr])

        [~, A_mf_vanilla.ITI] =...
            GenerateLatencyData_VanillaAlpha(params, A_mf_vanilla,...
            true, 'logn', 4);
        [ltom_vanilla{aa}(rr,:), htom_vanilla{aa}(rr,:)] =...
            block_dynamics_latency(A_mf_vanilla, twin, smoothfactor);

        [~, A_mf_rpe.ITI] =...
            GenerateLatencyData_RPEAlpha(params, A_mf_rpe,...
            true, 'logn', 4);
        [ltom_rpe{aa}(rr,:), htom_rpe{aa}(rr,:)] =...
            block_dynamics_latency(A_mf_rpe, twin, smoothfactor);
    end
end


%% Plot

figure
for aa = 1:3
    subplot(2, 3, aa); hold on
    plot(-twin:twin, mean(ltom_vanilla{aa}),'b')
    plot(-twin:twin, mean(htom_vanilla{aa}),'r')

    yl = ylim + [-0.015 0.05];
    fill([0 twin twin 0], [yl(1) yl(1) yl(2) yl(2)],...
        'k', facealpha=0.15, edgecolor='none')
    ylim(yl)

    xlim([-20 twin])

    if aa == 1
        title(['Vanilla model, alpha=' num2str(alphas(aa))])
    else
        title(['alpha=' num2str(alphas(aa))])
    end

    subplot(2, 3, aa+3); hold on
    plot(-twin:twin, mean(ltom_rpe{aa}),'b')
    plot(-twin:twin, mean(htom_rpe{aa}),'r')

    yl = ylim + [-0.015 0.015];
    fill([0 twin twin 0], [yl(1) yl(1) yl(2) yl(2)],...
        'k', facealpha=0.15, edgecolor='none')
    ylim(yl)

    xlim([-20 twin])
    if aa == 1
        title(['RPE*alpha model, alpha=' num2str(alphas(aa))])
    else
        title(['alpha=' num2str(alphas(aa))])
    end

end

set(gcf, position=[585 339 1096 606])
end