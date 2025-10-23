function [b, er] =...
    regress_modelRPE_vs_delay(A, RPE, nback, mixedblock)

if nargin<4
    mixedblock = 0; % use all blocks by default
end

if mixedblock==0
    usethese = A.hits==1;
else
    usethese = A.hits==1 & A.block==1;
end

A.reward_delay(A.reward_delay==100) = nan;
RD = log2(A.reward_delay);
RD(A.hits~=1) = nan;

% Create design matrix. This includes post-violation trials.
x = nan(length(RD), nback+2);
for row=1:length(RD)
    recent = RD(1:row);
    recent = recent(~isnan(recent));
    recent = recent(end:-1:max(1,end-nback));
    x(row, 1:length(recent)) = recent;
end
x(:,end) = ones(length(RD),1); % constant term.
X = x(usethese,:);
y = RPE(usethese);
[b, bint] = regress(y,X); % do regression.
er = bint(:,2)-bint(:,1); % calculate confidence intervals.
b = b(1:nback);
er = er(1:nback);

end