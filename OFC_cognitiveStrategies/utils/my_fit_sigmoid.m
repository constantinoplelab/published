function [finalParams, MSE, initCond, fitParams] = my_fit_sigmoid(x, y, iters)
%uses fmincon to fit sigmoidal to data specified by 'x' and 'y'.
%inputs: x = vector of x values.
%        y = vector of y values.
%        iter = iterations of random initializations. 10-20 is reasonable.

lb = [-5, 0, 0, -5];
ub = [5, 5, 30, 5];
nparams = length(lb);

% params(1): upper asymtote
% params(2): slope
% params(3): x offset

% mySigmoid =...
%     @(x, param)...
%     param(1)./(1 + exp(-param(2)*(x-param(3))));
mySigmoid =...
    @(x, param)...
    param(4) + (param(1)-param(4))./(1 + exp(-param(2)*(x-param(3))));

objFun = @(param) sum((y - mySigmoid(x, param)).^2);

[fitParams, initCond] = deal(nan(iters, nparams));
MSE = nan(iters, 1);

myOptions = optimset('Display', 'off',...
    'DiffMinChange', 0.0001, ...
    'MaxIter', 1000, 'GradObj', 'off', ...
    'Algorithm', 'interior-point', 'TolX', .00001);

for ii = 1:iters
    %randomly select starting parameters.
    initCond(ii,:) = rand(1, nparams).*(ub-lb) + lb;
    
    [fitParams(ii,:), MSE(ii), ~, ~, ~] = ...
        fmincon(objFun, initCond(ii,:), [], [], [], [],...
        lb, ub, [], myOptions);
end

[~, i] = min(MSE);
finalParams = fitParams(i,:);

end



