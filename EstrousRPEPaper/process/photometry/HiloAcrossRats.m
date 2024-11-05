function [high_rats, low_rats, T] = HiloAcrossRats(thirdrew_arg,...
    ratlist, Pstruct, Bstruct)  

high_rats = cell(length(ratlist), 1);
low_rats = cell(length(ratlist), 1);

for rat = 1:length(ratlist)

    ratname = ratlist{rat};
    disp(ratname)

    [high, low, T] = hilo_all_events({ratname},...
        Bstruct.(ratname), Pstruct.(ratname), thirdrew_arg);

    high_rats{rat} = high;
    low_rats{rat} = low;

end

end