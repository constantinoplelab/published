function [data_bc] = baselineCorrect(T, t1bl, t2bl, data)
    
    [~ , t1] = min(abs(T-t1bl));
    [~ , t2] = min(abs(T-t2bl));
    
    baseline = min(mean(data(:,t1:t2), 'omitnan'));
    data_bc = data - baseline;
end
