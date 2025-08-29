function b = regressVolume(data, trials, varargin)


if strcmpi('low', varargin)
    use = find(data.block == 3 & data.optout == 1 & ...
        data.catch == 1);
elseif strcmpi('high', varargin)
    use = find(data.block == 2 & data.optout == 1 & ...
        data.catch ==1);
else
    use = find(data.block==1 & data.optout==1 & ...
        data.catch==1);
end

if ~isempty(trials)
    use = intersect(use, trials);
end
    
%control
r = log2(data.reward(use));
x(:,1) = ones(length(use),1); 
x(:,2) = r;

y = data.wait_time(use);

%regression
b = regress(y, x);

