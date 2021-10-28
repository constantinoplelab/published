function ordering = time2peak(X,ops)
%gives sample ordering by time to peak

[nn,nt] = size(X);

%ordering, per alignment
peakrange_start = [1,nt];
peaktimes = zeros(1,nn); %time to peak

if isfield(ops,'useZScore')
    if ops.useZScore
        data_zscore_k = zeros(nn,nt); %for data
        for j = 1:nn
            datj = X(j,:);
            datmean = mean(datj);
            datvar = var(datj).^(0.5);
            data_zscore_k(j,:) = (datj-datmean)/(datvar);   
        end
        X = data_zscore_k;
    end
end


for j = 1:nn
    [~,pind] = max(X(j,peakrange_start(1):peakrange_start(2)));
    peaktimes(j) = pind;
end
[~,ordering] = sort(peaktimes);



