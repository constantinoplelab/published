function [b, tau, er, exp_t_1_coeff, stats, bestfit, bint] =...
    regress_latency_vs_rew(A, nback, plotarg,...
    z_arg, hit_arg, blk_arg, varargin)
% regress_latency_vs_rew - regresses trial initiation time against previous
%   rewards, a la Wang/Uchida 2013.
% INPUTS: 
%   A - Rat behavioral data. E.g., BestFitEarly.(RAT_ID).All.ratTrial
%   nback - number of previous rewards to include as regressors.
%   plotarg - whether to plot regression coefficients for each rat
%   z_arg - whether to z-score trial initiation times
%   hit_arg - whether to only include successfully completed trials (hits)
%   varargin - additional conditions for which trials to use (e.g., early
%   vs. late)
% OUTPUTS:
%   b - regression coefficients
%   tau - tau of exponential fit to regression coefficients
%   er - Error for each coefficient
%   exp_t_1_coeff - exponential fit
%   stats - regression statistics (from regress.m)
%   bestfit - full exponential parameters
%   bint - confidence interval for each coefficient

% Decide which data to use
if isempty(varargin)
    usethese = true(size(A.block));
else
    usethese = varargin{1};
end

if nargin<6
    blk_arg = 0; % use mixed block only by default
end

%z-score latencies
L = A.ITI(usethese);
if z_arg==1
    L(L>prctile(L, 99)) = nan;
    L = (L - mean(L, 'omitnan'))./std(L, 'omitnan');
end

% Pull rewards (non-hit trials are 0 if hit_arg = True)
if hit_arg==1
    R = log2(A.reward).*A.hits;
else
    R = log2(A.reward);
end

% Which block to use.
if blk_arg==0
    R(A.block~=1) = nan;
end

%Create design matrix. This includes post-violation trials.
x = nan(length(R), nback+1);
for j = 1:nback
    x(:,j) = [zeros(j,1); R(1:end-(j))];    
end
x(:,end) = ones(length(R),1); %constant term.

X = x(usethese,:); %mixed blocks only

[b, bint, ~, ~, stats] = regress(L,X); %do regression.
er = bint(:,2)-bint(:,1); %calculate confidence intervals.


% Find the first regression coefficient not significantly different from 
% zero (first one with 0 in 95% CI) gets set to zero, and all coefficients 
% including and after that
i = find(bint(:,1) < 0 & bint(:,2) > 0, 1);
b(i:end-1) = 0; 
er(i:end-1) = 0; 
bint(i:end-1,:) = 0;

%fit exponential
[~, bestfit, ~] = fit_exp_decay(1:nback, b(1:nback)', 10,...
    [-10 0], [0 10]); 
f = @(x, bestfit) (bestfit(1)*exp(-x./bestfit(2)));
newx = 1:.1:nback;
newy = f(newx, bestfit);

if all(b(1:end-1)==0)
    disp(['all betas are zero, the first beta to have a CI that crosses zero is t-' num2str(i)])
end

if b(1)~=0 && bestfit(2) == 0 %if the coefficient for t-1 is not zero and tau is zero, set tau to NaN
    tau = NaN; %set output tau to NaN
    bestfit(2) = NaN; %set tau exponential parameter to NaN
    exp_t_1_coeff = f(1, bestfit);
elseif all(b(1:end-1)==0) && bestfit(2)~=0 %if all coefficients are zero, set tau and the coefficient for t-1 to NaN
    tau = NaN;
    bestfit(2) = NaN; 
    exp_t_1_coeff = NaN;
else
    tau = bestfit(2);
    exp_t_1_coeff = f(1, bestfit); %where it crosses t-1 b/c when b(2)==0, A went to lower bound
end

%makeplot.
if plotarg
    figure;
    line([0 nback+1], [0 0], 'Color', [0 0 0]);
    plot(1:nback, b(1:nback), '.k', 'MarkerSize', 10); hold on
    errorbar(1:nback, b(1:nback), er(1:nback),...
        'Color', [0 0 0], 'LineStyle', 'none');
    plot(newx, newy, 'r');
    set(gca, 'TickDir', 'out'); box off;
    xlabel('N trials back');
    ylabel('Regression Coefficient');
    title(strcat(['Latency Tau = ', num2str(tau)]));
    subtitle(['Latency A = ' num2str(exp_t_1_coeff)])
end
end