function [DA,  AUC, pval, rcorr, T, ITI_legend, bins]...
    = ITI_effect(bstruct, pstruct, ITIbins_density,...
    event, binnum, window)

bdata = bstruct.(event);
pdata = pstruct.(event); 

% define bins
ITIs = bdata.ITI;
ITIs(ITIs>prctile(ITIs,99)) = NaN; %remove outliers
ITIs = [ITIs(2:end); NaN]; %shift initiation times to follow trial
bins = log2(ITIbins_density);
ITI_legend = cell(binnum, 1);
for ka = 1:length(bins)
    ITI_legend{ka} = num2str(2.^bins(ka)); %go from log2 back to linear scale
end

%initiate variables
DA = cell(1, binnum);
T = linspace(-5, 10, size(pstruct.('CPIn'), 2)-1);

[~ , tzero] = min(abs(T+0));
[~ , tafter] = min(abs(T-window));

%index trial types
these = find(bdata.Block==1 & bdata.TrialType==2 & bdata.PrevTrialType~=2); %  & bdata.Catch==0 & bdata.OptOut==0
data = pdata(these, 2:end);

%find ITIs
iti = log2(ITIs);
itis = iti(these);

y = nan(1,length(bins));
AUC = nan(1,length(bins));
for j = 2:length(y)
    ix = itis>bins(j-1) & itis<=bins(j);
    pdata = data(ix,:);
    da_mean = mean(pdata, 'omitnan');
    DA{1, j} = da_mean;
    if window > 0
        AUC(1,j) = trapz(T(:,tzero:tafter),...
            da_mean(:,tzero:tafter));
    elseif window < 0
        AUC(1,j) = trapz(T(:,tafter:tzero),...
            da_mean(:,tafter:tzero));
    end
end

%compute correlation
AUC_bytrial = NaN(size(data,1),1);
for trial = 1:size(data,1)
    if window > 0
        AUC_bytrial(trial, 1) = trapz(T(:,tzero:tafter),...
            data(trial,tzero:tafter));
    elseif window < 0
        AUC_bytrial(trial, 1) = trapz(T(:,tafter:tzero),...
            data(trial,tafter:tzero));
    end
end
[R,P] = corrcoef(itis(1:end-1), AUC_bytrial(1:end-1), 'rows','complete');
pval = P(2,1);
rcorr = R(2,1);

end
