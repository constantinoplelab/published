function [DA_sub,  DA_err, T, tstartbl, tzero, tafter, yint, slope] =...
    rew_effect_stages(bstruct, pstruct, block, Alignments,...
    rewards, Stages, window)

%initiate variables
DA_sub = cell(length(rewards), length(Stages));
DA_err = cell(length(rewards), length(Stages));
yint = NaN(1, length(Stages)); 
slope = NaN(1, length(Stages)); 

T = linspace(-5, 10, size(pstruct.(Alignments{1}), 2)-1);
[~ , tzero] = min(abs(T+0));
[~ , tafter] = min(abs(T-window));
[~ , tstartbl] = min(abs(T+0.05));

bdata = bstruct.(Alignments{1});
numgrp = NaN(length(Stages), 1);
stages = bdata.stage;

for s = 1:length(Stages)

    %find this stage
    thisstage_index = cellfun(@(x) logical(sum(strcmp(x, Stages{s}))), stages);

    %save number of sessions in this stage group
    numgrp(s) = length(unique(bdata.UniqueDay(thisstage_index)));

    Rewards = convertreward(bdata.Reward);
    y = nan(1, length(rewards)); %initialize AUC output

    for rew = 1:length(rewards)

        %subset trial types (no post violations)
        these = Rewards==rew & thisstage_index &...
            bdata.Block==block & bdata.PrevTrialType~=2;

        pdata = pstruct.(Alignments{1})(these, 2:end);

        da_mat = mean(pdata, 'omitnan');
        baseline = min(da_mat(:,tstartbl:tzero));
        data_bc = da_mat - baseline;

        this_sem = std(pdata, 'omitnan')./sqrt(size(pdata, 1));
        err = [data_bc-this_sem...
            fliplr(data_bc+this_sem)];
        err(isnan(err)) = 0;
        DA_sub{rew, s} = data_bc;
        DA_err{rew, s} = err;
        
        y(rew) = trapz(T(:,tzero:tafter), data_bc(:,tzero:tafter));
    end
    const = ones(length(y), 1);
    X = [const, (1:length(rewards))'];  %design matrix
    [beta,~,~,~,~] = regress(y', X); %regression model
    yint(1, s) = beta(1);
    slope(1, s) = beta(2);
end


end
