function [BetasPro, BetasDi] = DARegressRewardHistory_estrous(ratlist,...
    Pstruct, Bstruct, auc1, auc2, nback, Stages)
%%% Dopamine AUC as a function of late alpha model-predicted RPE
% Just mixed blocks
% Note: robustfit produces longer taus than regress

%Get DA AUC x reward history betas, split by stage for mixed blocks 
BetasPro = NaN(length(ratlist), nback+2); %for robust fit
BetasDi = NaN(length(ratlist), nback+2); %for robust fit
for rat = 1:length(ratlist)
    
    %identify rat
    ratname = ratlist{rat};
    disp(ratname)

    %load bdata and pdata
    bdata = Bstruct.(ratname).('CPIn'); %get offer cue behavioral data
    pstruct = Pstruct.(ratname);

    %%make reward history matrix
    %set up time vector and time points
    T = linspace(-5, 10, size(pstruct, 2)-1); 
    %set up AUC window
    [~, i1] = min(abs(T - auc1));
    [~, i2] = min(abs(T - auc2));
    dt = mean(diff(T)); %change in time
    %set up cell arrays to populate
    dates = unique(bdata.UniqueDay);
    [R_cell, AUC_cell] = deal(cell(length(dates), 1));
    for dd = 1:length(dates)

        thisday = ismember(bdata.UniqueDay, dates(dd));
        trialNos = bdata.TrialNumber(thisday);
        R_cell{dd} = nan(max(trialNos), 1);
        AUC_cell{dd} = nan(max(trialNos), 1);
        R_cell{dd}(trialNos) = bdata.Reward(thisday);

        %get AUC
        AUC_cell{dd}(trialNos) = sum(pstruct(thisday, i1:i2), 2)*dt;
    end

    R_vector = cell2mat(R_cell);
    R_log = log2(R_vector);
    AUC = cell2mat(AUC_cell);
    X = nan(length(R_log), nback+1);
    for ii = 1:nback+1
        X(:,ii) = [nan(ii-1, 1); R_log(1:end-(ii-1))];
    end
    %remove trials not recorded in bdata    
    AUC(isnan(X(:,1))) = []; 
    X(isnan(X(:,1)),:) = []; 

    for s = 1:length(Stages)
        %subset to stage and mixed blocks
        %find this stage
        thisstage_index = cellfun(@(x) logical(sum(strcmp(x, Stages{s}))),...
            bdata.stage);
        usethese = thisstage_index; % & bdata.TrialType~=2 & bdata.Block==1
        if s == 1
            BetasPro(rat, :) =...
                robustfit(X(usethese, :), AUC(usethese))'; %b1 is intercept
        elseif s == 2
            BetasDi(rat, :) =...
                robustfit(X(usethese, :), AUC(usethese))'; %b1 is intercept
        end
    end
end

end

