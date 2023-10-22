function [b, tau, er] =...
    regress_latency_vs_rew(A, nback, plotarg, z_arg, varargin)
%regress latency to cpoke against rewards, a la Wang/Uchida 2013.
%also fits exponential. cmc 12/2/20
%inputs: A struct. If you give it the S struct, will convert it to A struct.
%nback specifies number of previous trials to regress.
%cmc 11/03/21.
% avm 23 Nov 21 - edit to only use rewards and latencies from mixed blocks

if nargin == 2
    plotarg = false;
end

if isfield(A, 'pd') %if it's an S struct.
    [A, ~, ~] = parse_data_from_mysql(A);
end

if isempty(varargin)
    usetheseExtra = true(size(A.block));
else
    usetheseExtra = varargin{1};
end

usethese = A.block == 1 & usetheseExtra;

%z-score latencies.
L = A.ITI(usethese);
if z_arg==1
    L(L>prctile(L, 99)) = nan;
    L = (L - mean(L, 'omitnan'))./std(L, 'omitnan');
end

R = log2(A.reward).*A.hits;
R(A.block~=1) = nan;

%Create design matrix. This includes post-violation trials.
x = nan(length(R), nback+1);
for j = 1:nback
    x(:,j) = [nan(j,1); R(1:end-(j))];    
end
x(:,end) = ones(length(R),1); %constant term.

X = x(usethese,:); %mixed blocks only

[b, bint] = regress(L,X); %do regression.
er = bint(:,2)-bint(:,1); %calculate confidence intervals.

% %this is just to get the fmincon code to fit better. 
% % values were really small after z-scoring.
% b = b.*10; 

% Find the first non-significant regression coefficient (first one with
% 0 in 95% CI) then zero out all coefficients including and after that
i = find(bint(:,1) < 0 & bint(:,2) > 0, 1);
b(i:end-1) = 0;
er(i:end-1) = 0;

[~, bestfit, ~] = fit_exp_decay(1:nback, b(1:nback)', 10,...
    [-10 0], [0 10]); %fit exponential
f = @(x, bestfit) (bestfit(1)*exp(-x./bestfit(2)));
newx = 1:.1:nback;
newy = f(newx, bestfit);

if all(b(2:end-1)==0)
    tau = 0;
else
    tau = bestfit(2);
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
end
end