function [DA, DA_err, T, AUC] = rew_effect(bstruct,...
    pstruct, block, rewards, Alignments, port)

DA = cell(length(rewards), length(Alignments));
DA_err = cell(length(rewards), length(Alignments));
AUC = nan(length(Alignments), length(rewards)); %initialize AUC output

for a = 1:length(Alignments)
    
    T = linspace(-5, 10, size(pstruct.(Alignments{a}), 2)-1);

    bdata = bstruct.(Alignments{a});
    Rewards = convertreward(bdata.Reward);

    for r = 1:length(rewards)
        
        rew = rewards(r);

        %subset trial types
        if strcmp(Alignments{a}, 'CPIn')
            if port==0
                these = Rewards==rew & bdata.Block==block &...
                    bdata.PrevTrialType~=2;
            else
                these = Rewards==rew & bdata.Block==block &...
                    bdata.PrevTrialType~=2 & bdata.RewardPort==port;
            end
        elseif strcmp(Alignments{a}, 'OptOut')
            if port==0
                these = Rewards==rew & bdata.OptOut==1 & bdata.Block==block...
                    & bdata.PrevTrialType~=2;
            else
                these = Rewards==rew & bdata.OptOut==1 & bdata.Block==block...
                    & bdata.PrevTrialType~=2 & bdata.RewardPort==port;
            end
        else
            if port==0
                these = Rewards==rew & ...
                    bdata.Block==block & bdata.PrevTrialType~=2;
            else
                these = Rewards==rew & ...
                    bdata.Block==block & bdata.PrevTrialType~=2 ...
                    & bdata.RewardPort==port;
            end
        end
        pdata = pstruct.(Alignments{a})(these, 2:end);                

        da_mean = mean(pdata, 'omitnan');
        standarderr = sem(pdata);
        DA{rew, a} = da_mean;  
        DA_err{rew, a} = standarderr;
        
    end
end


