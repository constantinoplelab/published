function [DA, AUC, yint, slope] = delay_effect_stages(bstruct, ...
    pstruct, block, event, delaybins, window, Stages)
%block = 1, 2, or 3
%delays_rats = empirically determined delay bins

yint = NaN(1, 2); 
slope = NaN(1, 2); 
bins = delaybins;

%initiate variables
DA = cell(length(Stages), length(bins));
AUC = NaN(length(Stages),length(bins));    
T = linspace(-5, 10, size(pstruct.(event), 2)-1);
[~ , tzero] = min(abs(T+0));
[~ , tbefore] = min(abs(T+0.05));
[~ , tafter] = min(abs(T-window));
bdata = bstruct.(event);

%only use staged data, also create stages vector
stages = bdata.stage;
numgrp = NaN(length(Stages), 1);
for s = 1:length(Stages)
    %find stage
    thisstage_index = cellfun(@(x) logical(sum(strcmp(x, Stages{s}))),...
        stages);
    %save number of sessions in this stage group
    numgrp(s) = length(unique(bdata.UniqueDay(thisstage_index)));
    %index trial types
    these = find(bdata.Block==block & bdata.Catch==0 &...
        bdata.OptOut==0 & thisstage_index);  
    data = pstruct.(event)(these, 2:end);
    %find delay
    delay = bdata.RewardDelay(these);
    delay_legend = cell(length(bins), 1);
    for ka = 1:length(bins)
        delay_legend{ka} = num2str(bins(ka));
    end
    d = delay;

    for j = 2:length(AUC)
        ix = find(d>bins(j-1) & d<=bins(j));
        da_mat = mean(data(ix,:), 'omitnan'); 
        da_mat_bl = da_mat - min(mean(data(ix,tbefore:tzero), 'omitnan'));
        if length(da_mat)>1
            DA{s, j} = da_mat; %save uncorrected version
            if window<0
                AUC(s, j) = trapz(T(:,tafter:tzero), da_mat_bl(:,tafter:tzero)); %AUC
            elseif window>0
                AUC(s, j) = trapz(T(:,tzero:tafter), da_mat_bl(:,tzero:tafter)); %AUC
            end
        else
           DA{s, j} = NaN(1, length(T));
           AUC(s, j) = NaN;
        end
    end

    %get slope and intercept for AUCs by delay bin
    x = bins(1:end-1);
    const = ones(length(AUC(s, 2:end)),1);   
    X = [const, x];  %design matrix
    [beta,~,~,~,~] = regress(AUC(s, 2:end)', X); %regression model
    yint(1, s) = beta(1); 
    slope(1, s) = beta(2);    

end


end
