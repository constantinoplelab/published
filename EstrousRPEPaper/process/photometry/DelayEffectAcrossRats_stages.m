function [da_rats, AUC_rats, yint_rats, slp_rats] =...
    DelayEffectAcrossRats_stages(ratlist, block, delaybins,...
    window, Pstruct, Bstruct, event, Stages)
%Sorts response to SideOff by delay length
%input
    %block = behavior block to use
    %numbins = number of delay bins

da_rats = cell(length(ratlist), 1);
AUC_rats = cell(length(ratlist), 1);
yint_rats = cell(length(ratlist), 1);
slp_rats = cell(length(ratlist), 1);

for rat = 1:length(ratlist)
    ratname = ratlist{rat};
    disp(ratname)
    [DA, AUC, yint, slope] =...
        delay_effect_stages(Bstruct.(ratname), Pstruct.(ratname), block,...
        event, delaybins, window, Stages);
    da_rats{rat} = DA;
    AUC_rats{rat} = AUC;
    yint_rats{rat} = yint;
    slp_rats{rat} = slope;
    
end

end
