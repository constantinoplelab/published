function [initcond, bestfit, mse] = fit_exp_decay(x, y, iters, varargin)
%uses fmincon to fit exponential to data specified by 'x' and 'y'.
%inputs: x = vector of x values.
%        y = vector of y values.
%        iter = iterations of random initializations. 10-20 is reasonable.
%cmc 11/3/21
%adding upper and lower bounds as varargin cmc 7/28/23

if isempty(varargin)
    lower_bound = [-50 -50];
    upper_bound = [50 50];
else
    lower_bound = varargin{1};
    upper_bound = varargin{2};
end

data.x = x;
data.y = y;
nparams = 2;

initcond = rand(iters, nparams).*(upper_bound-lower_bound) + lower_bound;
f_out = initcond;
mse = nan(iters,1);

for mx = 1:iters
   
    [f_out(mx,:), mse(mx), exitflag, ~, ~, ~] =...
        fitfun(data, initcond(mx,:), lower_bound, upper_bound);
    
end

[~, i] = min(mse);
bestfit = f_out(i,:);


end

function [x_fmincon, f, exitflag, output, grad, hessianmat] = fitfun(data, ...
    x_init, lower_bound, upper_bound)

obj_fun = @(x)fit_hyper_fun(x, data);

[x_fmincon, f, exitflag, output, ~, grad, hessianmat] = ...
    fmincon(obj_fun, x_init, [], [], [], [], lower_bound, upper_bound, [], ...
    optimset('Display', 'off', ...
    'DiffMinChange', 0.0001, ...
    'MaxIter', 500, 'GradObj', 'off', ...
    'Algorithm', 'interior-point', 'TolX', .00001));
end


function [minimized]=fit_hyper_fun(param, data)

f = @(x) (param(1)*exp(-x./param(2)));

y = f(data.x);
minimized = nanmean((data.y-y).^2); %minimize mean squared error

end
