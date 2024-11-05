function [DAbHilo, stageffect, pval] =...
    HiloAcrossRats_stages(window, thirdrew_arg, ratlist,...
    Pstruct, Bstruct, Stages)
%%for plotting response aligned to event in high versus low blocks across stages, 
    %averaged over rats

DAbHilo = [];

for rat = 1:length(ratlist)

    ratname = ratlist{rat};
    disp(ratname)

    [high, low, ~, ~, deltas] = hiloRew_stages(Bstruct.(ratname),...
        Pstruct.(ratname), window, thirdrew_arg, Stages);

    DAbHilo.high{rat} = high;
    DAbHilo.low{rat} = low;
    DAbHilo.deltas{rat} = deltas;

end

%get stage effect (proestrus - diestrus)
if length(ratlist) > 1
    deltas_rats = NaN(length(ratlist), 2);
    for rat = 1:length(ratlist)
        deltas_rats(rat, :) = DAbHilo.deltas{rat}';
    end
    pval = signrank(deltas_rats(:, 1), deltas_rats(:, 2));
    stageffect = deltas_rats(:, 1) - deltas_rats(:, 2);
end

end
