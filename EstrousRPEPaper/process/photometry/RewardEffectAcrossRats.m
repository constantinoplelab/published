function [da_rats, da_err_rats, T] = RewardEffectAcrossRats(block,...
    Alignments, port, ratlist, Pstruct, Bstruct) %rmv_postvios
%use varargin to input rat list

if block==0
    rewards = [1 2 3 4 5];
elseif block==1
    rewards = [1 2 3 4 5];
elseif block==2
    rewards = [3 4 5];
elseif block==3
    rewards = [1 2 3];
end

%initiate variables
da_rats = cell(length(ratlist), 1);
da_err_rats = cell(length(ratlist), 1);
AUC_y_rats = cell(length(ratlist), 1);

for rat = 1:length(ratlist)    
    
    ratname = ratlist{rat};

    disp(ratname)

    %Load data and run function
    [DA, DA_err, T, AUC] = rew_effect(Bstruct.(ratname),...
        Pstruct.(ratname), block, rewards, Alignments, port);

    da_rats{rat} = DA;
    da_err_rats{rat} = DA_err;
    AUC_y_rats{rat} = AUC;  
        
end

end
