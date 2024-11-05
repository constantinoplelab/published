function [AUC_overgroup, DA] = rew_effect_byblock_stages(bstruct,...
    pstruct, window, Stages) 

AUC_overgroup = [];
DA = cell(2,1);
blocks = 1:3;
allrewards = [4 8 16 32 64];

%getdata
bdata = bstruct.('CPIn');
stages = bdata.stage;
AUC = NaN(length(blocks)*2, length(allrewards)); %initialize AUC output for five rewards
T = linspace(-5, 10, size(pstruct.('CPIn'), 2)-1);
numsess_stages = NaN(length(Stages) , 1);

for s = 1:length(Stages)
    DA_mat = cell(3, 5);

    %find this stage
    thisstage_index = cellfun(@(x) logical(sum(strcmp(x, Stages{s}))), stages);
    numsess_stages(s) = length(unique(bdata.UniqueDay(thisstage_index)));

    for bl = 1:length(blocks)
    
        if bl==1
            rewards = [1 2 3 4 5];
        elseif bl==2
            rewards = [3 4 5];
        elseif bl==3
            rewards = [1 2 3];
        end
            
        [~ , tzero] = min(abs(T+0));
        [~ , tstartbl] = min(abs(T+0.05));
        [~ , tafter] = min(abs(T-window));
    
        Rewards = convertreward(bdata.Reward);
    
        for r = 1:length(rewards)
            
            rew = rewards(r);
            
            %subset trial types
            these = Rewards==rew & bdata.Block==bl & thisstage_index...
                & bdata.PrevTrialType~=2; %exclude post violations
            pdata = pstruct.('CPIn')(these, 2:end); 
    
            da_mat = mean(pdata, 'omitnan');
            baseline = min(da_mat(:,tstartbl:tzero));
            data_bc = da_mat - baseline;
            DA_mat{bl, r} = data_bc;
            AUC(length(blocks)*(s-1)+bl, rew) = trapz(T(:,tzero:tafter),...
                data_bc(:,tzero:tafter));
        end
    end
    DA{s} = DA_mat;
end

%Normalize to mixed 64 ul
AUC_normalized = (AUC-min(AUC(:)))./(max(AUC(:)) - min(AUC(:)));
AUC_overgroup.(Stages{1}).AUC = AUC_normalized(1:3,:);
AUC_overgroup.(Stages{2}).AUC = AUC_normalized(4:6,:);

end
