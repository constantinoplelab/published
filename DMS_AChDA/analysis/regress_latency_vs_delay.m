function [b, er] = regress_latency_vs_delay(A, nback, z_arg, mixedblock)

if nargin<4
    mixedblock = 0; % use all blocks by default
end

% filter outliers using all trials, not just the subset given
L = A.ITI;
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

posthit = [0; A.hits(1:end-1)];
if mixedblock==0
    usethese = posthit==1;
else
    usethese = posthit==1 & A.block==1;
end

A.reward_delay(A.reward_delay==100) = nan;
RD = log2(A.reward_delay);
RD(A.hits~=1) = nan;

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
y = L(usethese);
[b, bint] = regress(y,X); %do regression.
er = bint(:,2)-bint(:,1); %calculate confidence intervals.
b = b(1:nback);
er = er(1:nback);

end