function [AUC_overgroup, DA, DA_normalized_mixed] = rew_effect_byblock_stages(bstruct,...
    pstruct, window, Stages) 

%set variables
AUC_overgroup = [];
DA = cell(length(Stages),1);
DA_normalized_mixed = cell(length(Stages),1);
blocks = 1:3;
allrewards = [4 8 16 32 64];

%getdata
bdata = bstruct.('CPIn');
stages = bdata.stage;
AUC = NaN(length(blocks)*length(Stages), length(allrewards)); %initialize AUC output for five rewards
if size(pstruct, 2)>1
    T = linspace(-5, 10, size(pstruct, 2)-1);
else
    T = linspace(-5, 10, size(pstruct.('CPIn'), 2)-1);   
end
numsess_stages = NaN(length(Stages) , 1);

for s = 1:length(Stages)
    DA_mat = cell(length(blocks), length(allrewards));

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
            if size(pstruct,1) > 1
                pdata = pstruct(these, 2:end); 
            else
                pdata = pstruct.('CPIn')(these, 2:end);                 
            end
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

%Min-max normalize AUC
AUC_normalized = (AUC-min(AUC(:)))./(max(AUC(:)) - min(AUC(:)));
if length(Stages)==2
    AUC_overgroup.(Stages{1}).AUC_norm = AUC_normalized(1:3,:);
    AUC_overgroup.(Stages{2}).AUC_norm = AUC_normalized(4:6,:);
    AUC_overgroup.(Stages{1}).AUC = AUC(1:3,:);
    AUC_overgroup.(Stages{2}).AUC = AUC(4:6,:); 
elseif length(Stages)==4
    AUC_overgroup.(Stages{1}).AUC_norm = AUC_normalized(1:3,:);
    AUC_overgroup.(Stages{2}).AUC_norm = AUC_normalized(4:6,:);
    AUC_overgroup.(Stages{3}).AUC_norm = AUC_normalized(7:9,:);
    AUC_overgroup.(Stages{4}).AUC_norm = AUC_normalized(10:12,:);
    AUC_overgroup.(Stages{1}).AUC = AUC(1:3,:);
    AUC_overgroup.(Stages{2}).AUC = AUC(4:6,:); 
    AUC_overgroup.(Stages{3}).AUC = AUC(7:9,:);
    AUC_overgroup.(Stages{4}).AUC = AUC(10:12,:);    
end

end
