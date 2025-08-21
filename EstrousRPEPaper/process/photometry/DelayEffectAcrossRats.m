function [da_rats, T, bins, delay_legend] = DelayEffectAcrossRats(ratlist,...
    block, delaybins, event, Bstruct, Pstruct)
%input
    %block = behavior block to use
    %numbins = number of delay bins

da_rats = cell(length(ratlist), 1);

for rat = 1:length(ratlist)    
    
    ratname = ratlist{rat};

    disp(ratname)

    [DA,  T, ~, delay_legend, bins] =...
        delay_effect_cg(Bstruct.(ratname), Pstruct.(ratname),...
        event, block, delaybins);

    da_rats{rat} = DA;

end

end
