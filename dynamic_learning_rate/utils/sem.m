function s = sem(data, nan_arg)    
    if nargin<2
        nan_arg = 'omitnan';
    end
    s = std(data, nan_arg)./ sqrt(size(data, 1));
end