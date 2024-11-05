function [b, tau, er, exp_t_1_coeff, stats, bestfit, bint] =...
    regress_latency_vs_rew(A, nback, plotarg,...
    z_arg, hit_arg, blk_arg, varargin)
%regress latency to cpoke against rewards, a la Wang/Uchida 2013.
%also fits exponential. cmc 12/2/20
%inputs: A struct. If you give it the S struct, will convert it to A struct.
%nback specifies number of previous trials to regress.
%cmc 11/03/21.
% avm 23 Nov 21 - edit to only use rewards and latencies from mixed blocks

if isfield(A, 'pd') %if it's an S struct.
    [A, ~, ~] = parse_data_from_mysql(A);
end

if isempty(varargin)
    usethese = true(size(A.block));
else
    usethese = varargin{1};
end

if nargin<6
    blk_arg = 0; % use mixed block only by default
end

%z-score latencies.
L = A.ITI(usethese);
if z_arg==1
    lb = prctile(L, 1);
    ub = prctile(L, 99);
    L(L < lb | L > ub) = nan;
    L = (L - mean(L, 'omitnan'))./std(L, 'omitnan');
end

if hit_arg==1
    R = log2(A.reward).*A.hits;
else
    R = log2(A.reward);
end

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

% %this is just to get the fmincon code to fit better. 
% % values were really small after z-scoring.
% b = b.*10; 

% Find the first regression coefficient not significantly different from zero (first one with
% % 0 in 95% CI) gets set to zero, and all coefficients including and after that
i = find(bint(:,1) < 0 & bint(:,2) > 0, 1);
b(i:end-1) = 0; 
er(i:end-1) = 0; 
bint(i:end-1,:) = 0;

[~, bestfit, ~] = fit_exp_decay(1:nback, b(1:nback)', 10,...
    [-10 0], [0 10]); %fit exponential
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