function [b, tau, er, exp_t_1_coeff, stats, bestfit, bint, rint] =...
    regress_latency_vs_rew(A, nback, z_arg, hit_arg, blk_arg, varargin)

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

L = A.ITI(usethese);

% filter outliers using all trials, not just the subset given
thresh = [1, 99];
cutoff = [prctile(A.ITI,min(thresh)), prctile(A.ITI,max(thresh))];
fprintf('ITI filtered at %.1f and %.1f percent: %.2f, %.2f\n', ...
    min(thresh), max(thresh), min(cutoff), max(cutoff))
L(L<min(cutoff) | L>max(cutoff)) = nan;
% de-trend
disp('De-trending..')
x = [1:length(L)]';
X = [x, ones(length(x),1)];
b = regress(L, X);
newy = b(1)*x + b(2);
L = L - newy;

if z_arg==1
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

X = x(usethese,:);

[b, bint, ~, rint, stats] = regress(L,X); %do regression.
er = (bint(:,2)-bint(:,1))/2; %calculate confidence intervals.


% % Find the first non-significant regression coefficient (first one with
% % 0 in 95% CI) then zero out all coefficients including and after that
% i = find(bint(:,1) < 0 & bint(:,2) > 0, 1);
% b(i:end-1) = 0;
% er(i:end-1) = 0;


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

end