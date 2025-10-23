function plot_fig4g(datadir)
% Plot time to opt-out regressed against reward delay and reward volume in
% the past 15 trials during mixed blocks, averaged acros rats. Ten trials
% are shown for visualization purposes. N = 16 rats.
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

load(fullfile(datadir, 'data-published/ratlist.mat'), 'ratList');
nback = 15;

figure;
set(gcf, units='inches', position=[7,4,4.5,2])
t = tiledlayout(1,2);

% opt-out time regressed against reward volumes
nexttile(1); hold on
b = [];
for rr=1:length(ratList)
    rat = ratList{rr};
    fprintf('%i of %i: %s\n', rr, length(ratList), rat)

    [coeff, er] = wt_regress_rew(datadir, rat, nback);

    if sum(isnan(er))==0 && sum(er==0)<nback
        b = [b, coeff];
    end
end

b = b';
plotPretty(1:10, b(:,1:10), '#793CAB')
xlim([0.5 10.5])
ylim([-0.05 0.25])
yticks(0:0.1:0.2)
yline(0, 'k--')
set(gca, xtick=[1:1:10], xticklabels=[0:1:9])
[~,pval] = ttest(b);
sig = find(pval<0.05);
if ~isempty(sig)
    arrayfun(@(x) fprintf('%i trial back: p=%.5f\n', x-1, pval(x)),...
        sig)
end

subtitle('vs reward volume');
set(gca, fontsize=8)

% opt-out time regressed against reward delays
nexttile(2); hold on
b = [];
for rr=1:length(ratList)
    rat = ratList{rr};
    fprintf('%i of %i: %s\n', rr, length(ratList), rat)

    [coeff, er] = wt_regress_delay(datadir, rat, nback);

    if sum(isnan(er))==0 && sum(er==0)<nback
        b = [b, coeff];
    end
end

b = b';

plotPretty(1:10, b(:,1:10), 'k')
xlim([0.5 10.5])
ylim([-0.05 0.25])
yticks(0:0.1:0.2)
yline(0, 'k--')
set(gca, 'xtick', 1:10)
[~,pval] = ttest(b);
disp(pval)
sig = find(pval<0.05);
if ~isempty(sig)
    arrayfun(@(x) fprintf('%i trial back: p=%.5f\n', x, pval(x)),...
        sig)
end

subtitle('vs reward delay');
set(gca, fontsize=8)

xlabel(t, 'Trials back', fontsize=8)
ylabel(t, 'Regression weight', fontsize=8)



end

function [b, er] = wt_regress_rew(datadir, rat, nback)
% regress wait time in mixed blocks against reward offers

folder = fullfile(datadir, 'data-published/A_structs');
file = strcat('ratTrial_', rat);
load(fullfile(folder, file), 'A');

usethese = A.optout==1 & A.block==1;

% z-score wait times using all opt-out trials
WT = A.wait_time;
WT(A.optout==0) = nan;
WT = (WT - mean(WT, 'omitnan'))./(std(WT, 'omitnan'));

R = log2(A.reward);
R(A.block~=1) = nan;

%Create design matrix. This includes post-violation trials.
x = nan(length(R), nback+1);
for j = 1:nback
    x(:,j) = [zeros(j-1,1); R(1:end-(j-1))];
end
x(:,end) = ones(length(R),1); %constant term.
X = x(usethese,:);
y = WT(usethese);
[b, bint] = regress(y,X); %do regression.
er = bint(:,2)-bint(:,1); %calculate confidence intervals.
b = b(1:nback);
er = er(1:nback);
end

function [b, er] = wt_regress_delay(datadir, rat, nback)
% regress wait time in mixed blocks against reward delays

folder = fullfile(datadir, 'data-published/A_structs');
file = strcat('ratTrial_', rat);
load(fullfile(folder, file), 'A');

posthit = [0; A.hits(1:end-1)];
usethese = A.optout==1 & A.block==1 & posthit==1;

% z-score wait times using all opt-out trials
WT = A.wait_time;
WT(A.optout==0) = nan;
WT = (WT - mean(WT, 'omitnan'))./(std(WT, 'omitnan'));

A.reward_delay(A.reward_delay==100) = nan;
RD = log2(A.reward_delay);
RD(A.block~=1 | A.hits~=1) = nan;

%Create design matrix
x = nan(length(RD), nback+1);
for row=2:length(RD)
    recent = RD(1:row-1);
    recent = recent(~isnan(recent));
    recent = recent(end:-1:max(1,end-nback));
    x(row, 1:length(recent)) = recent;
end
x(:,end) = ones(length(RD),1); %constant term.
X = x(usethese,:);
y = WT(usethese);
[b, bint] = regress(y,X); %do regression.
er = bint(:,2)-bint(:,1); %calculate confidence intervals.
b = b(1:nback);
er = er(1:nback);
end






