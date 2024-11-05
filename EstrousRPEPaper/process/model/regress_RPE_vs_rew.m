function [b, er, stats] = regress_RPE_vs_rew(A, nback)
%regress model-predicted RPEs against rewards
%inputs: A struct. If you give it the S struct, will convert it to A struct.
%nback specifies number of previous trials to regress.

if isfield(A, 'pd') %if it's an S struct.
    [A, ~, ~] = parse_data_from_mysql(A);
end

%get variables
R = log2(A.reward);
RPEs = A.RPE;

%Create design matrix. This includes the current reward.
X = nan(length(R), nback+1);
for j = 1:nback+1
    X(:,j) = [zeros(j-1,1); R(1:end-(j-1))];    
end
X(:,end) = ones(length(R),1); %constant term.

[b, bint, ~, ~, stats] = regress(RPEs,X); %do regression.
er = bint(:,2)-bint(:,1); %calculate confidence intervals.

% %this is just to get the fmincon code to fit better. 
% % values were really small after z-scoring.
% b = b.*10; 

% Find the first regression coefficient not significantly different from zero (first one with
% % 0 in 95% CI) gets set to zero, and all coefficients including and after that
i = find(bint(:,1) < 0 & bint(:,2) > 0, 1);
b(i:end-1) = 0; 
er(i:end-1) = 0; 

if all(b(1:end-1)==0)
    disp(['all betas are zero, the first beta to have a CI that crosses zero is t-' num2str(i)])
end

end