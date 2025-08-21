function [high_rats, low_rats, T] = HiloAcrossRats(thirdrew_arg,...
    ratlist, Pstruct, Bstruct)  

high_rats = cell(length(ratlist), 1);
low_rats = cell(length(ratlist), 1);

for rat = 1:length(ratlist)

    ratname = ratlist{rat};
    disp(ratname)

    [high, low] = hilo_one_event(Bstruct.(ratname), ...
        Pstruct.(ratname), thirdrew_arg, 'CPIn');

    high_rats{rat} = high;
    low_rats{rat} = low;

end

end