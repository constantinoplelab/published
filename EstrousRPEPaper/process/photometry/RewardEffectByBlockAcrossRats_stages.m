function AUC_byrat =...
    RewardEffectByBlockAcrossRats_stages(window, Stages, ratlist,...
    Pstruct, Bstruct)

%initiate variables
AUC_byrat = cell(1,length(ratlist));
for rat = 1:length(ratlist)
    ratname = ratlist{rat};
    disp(ratname)
    AUC_byrat{rat} =...
        rew_effect_byblock_stages(Bstruct.(ratname),...
        Pstruct.(ratname), window, Stages);
end

end


