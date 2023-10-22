function FigureS6(datadir, codedir)
%FigureS6 -  Satiety effects for trial initiation time are modest and do 
% not qualitatively affect results
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

%% Process data
a = load([datadir filesep 'ratList.mat']);
ratList = a.ratList; % list of rats

% block transition dynamic parameters
twin = 20; % number of trials around each block to pull
smoothfactor = 10; % size of smoothing window

% Trial initiation time session parameters
example_rat = 'L027'; % Wwhich rat to use
i = 3; % Which session to pull

% Pre-allocate matricies
[ltomNoDT, htomNoDT, mtolNoDT, mtohNoDT] =...
    deal(nan(length(ratList), 2*twin+1));

% block transition without detrending over rats
for rr = 1:length(ratList)
    fprintf('%d out of %d\n', rr, length(ratList));

    % Load Data
    fname = strcat(['ratTrial_', ratList{rr}, '.mat']);
    A = load([datadir 'A_Structs_Final' filesep fname]);
    A = A.A;

    if strcmp(ratList{rr}, example_rat)
        firstTrial = [find(A.trial_num==1); length(A.trial_num)+1];

        usethese = firstTrial(i):firstTrial(i+1)-1;

        trialNum = A.trial_num(usethese);
        L = A.ITI(usethese);
        L(1) = nan;
        L(L>prctile(L, 99)) = nan;

        betas = robustfit(trialNum, L);
    end

    % no detrending trial initiation time
    [ltomNoDT(rr,:), htomNoDT(rr,:),...
        mtolNoDT(rr,:), mtohNoDT(rr,:)] =...
        block_dynamics_latency(A, twin, smoothfactor, false);
end


%% Plot

figure
clear l

subplot(1, 3, 1); hold on
plot(trialNum, L, 'k.')
refline(betas(2), betas(1))
ylim([0 2])
xlim([0 410])

xlabel('Trial Number')
ylabel('Trial initiation time (s)')

title(['y = ' num2str(betas(2)) '*x + ' num2str(betas(1))])

subplot(1, 3, 2); hold on
fill([0 20 20 0], [-0.2 -0.2 0.2 0.2], 'k',...
    facealpha=0.1, linestyle='none')

l(1) = shadedErrorBar(-twin:twin, mean(ltomNoDT),...
    std(ltomNoDT)./sqrt(length(ratList)), Lineprops={'b'});
l(2) = shadedErrorBar(-twin:twin, mean(htomNoDT),...
    std(htomNoDT)./sqrt(length(ratList)), Lineprops={'r'});

ylim([-0.2 0.2])

ylabel('Trial initiation time (z-score)')

subplot(1, 3, 3); hold on

fill([-20 0 0 -20], [-0.2 -0.2 0.2 0.2], 'k',...
    facealpha=0.1, linestyle='none')

l(3) = shadedErrorBar(-twin:twin, mean(mtolNoDT),...
    std(mtolNoDT)./sqrt(length(ratList)), Lineprops={'b'});
l(4) = shadedErrorBar(-twin:twin, mean(mtohNoDT),...
    std(mtohNoDT)./sqrt(length(ratList)), Lineprops={'r'});

ylim([-0.2 0.2])

arrayfun(@(line) arrayfun(@(e) set(e, LineStyle='none'), line.edge), l)
set(gcf, position=[290 527 1390 420])

end