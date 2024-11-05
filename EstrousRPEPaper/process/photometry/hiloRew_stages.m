function [high, low, err_high, err_low, deltas, T, tstartbl, tendbl, tzero]...
    = hiloRew_stages(bstruct, pstruct, window, thirdrew_arg,...
        Stages) 
%%plot response aligned to offer cue in high versus low blocks, separated by
%%input:
    %ratlist = which rats to run through on at a time
    %pstruct = aligned photometry data by rat
    %bstruct = aligned behavioral data by rat    
    
high = cell(length(Stages), 1);
low = cell(length(Stages), 1);
err_high = cell(length(Stages), 1);
err_low = cell(length(Stages), 1);
deltas = NaN(length(Stages), 1);
block = {'mixed'; 'high'; 'low'};
bdata = bstruct.('CPIn');   
stages = bdata.stage;  
Rewards = convertreward(bdata.Reward);

T = linspace(-5, 10, size(pstruct.('CPIn'), 2)-1);
[~ , tzero] = min(abs(T));
[~ , tendbl] = min(abs(T+0));
[~ , tstartbl] = min(abs(T+0.05));
[~ , tafter] = min(abs(T-window));
numsess_stages = NaN(length(Stages) , 1);
for s = 1:length(Stages) 

    %find this stage
    thisstage_index = cellfun(@(x) logical(sum(strcmp(x, Stages{s}))), stages);
    numsess_stages(s) = length(unique(bdata.UniqueDay(thisstage_index)));
        
    bhi = NaN(1, size(T, 2));
    blo = NaN(1, size(T, 2));

    AUC = NaN(2, 1);

    for bl = 2:length(block)

        if thirdrew_arg==1
            these = bdata.Block==bl & Rewards==3 ...
                & thisstage_index & bdata.PrevTrialType~=2;
        else
            these = bdata.Block==bl ...
                & thisstage_index & bdata.PrevTrialType~=2;
        end

        %get pstruct data with these criteria
        pdata = pstruct.('CPIn')(these, 2:end);    

        da_mat = mean(pdata, 'omitnan');
        baseline = min(da_mat(:,tstartbl:tendbl));
        y = da_mat - baseline;

        sem = std(pdata, 'omitnan')./sqrt(size(pdata, 1));
        err = [y-sem...
            fliplr(y+sem)];
        err(isnan(err)) = 0;
        if bl == 2
            bhi(1, :) = y;
            hierr = err;
            AUC(1, 1) = trapz(T(:,tzero:tafter), y(:,tzero:tafter));
        elseif bl == 3
            blo(1, :) = y;
            loerr = err;
            AUC(2, 1) = trapz(T(:,tzero:tafter), y(:,tzero:tafter));
        end

    end
 
    high{s, 1} = bhi;
    low{s, 1} = blo;
    err_high{s ,1} = hierr;
    err_low{s, 1} = loerr;
    deltas(s, 1) = AUC(1, 1) - AUC(2, 1); %high - low

end       


end

