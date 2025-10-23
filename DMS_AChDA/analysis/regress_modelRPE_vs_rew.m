function [b] =...
    regress_modelRPE_vs_rew(A, RPE, nback, hit_arg, allblocks)

if nargin<5
    allblocks = 0; % use mixed block only by default
end

if hit_arg==1
    R = log2(A.reward).*A.hits;
else
    R = log2(A.reward);
end

if allblocks==0
    usethese = A.block==1;
    R(~usethese) = nan;
end

% Create design matrix. This includes post-violation trials.
x = nan(length(R), nback+2);
for j = 1:nback+1
    x(:,j) = [zeros(j-1,1); R(1:end-j+1)];    
end
x(:,end) = ones(length(R),1); % constant term.

X = x(usethese,:); % mixed blocks only
RPE(~usethese) = [];

[b, ~,~,~,~] = regress(RPE,X); %do regression.


end