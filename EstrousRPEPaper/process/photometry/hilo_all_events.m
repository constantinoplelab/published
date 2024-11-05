function [high, low, T, bhi, blo] = hilo_all_events(ratList, bstruct, pstruct,...
    thirdrew_arg)
%%plot response aligned to event in high versus low blocks
%%input:
%ratlist = which rats to run through on at a time
%pstruct = aligned photometry data by rat
%bstruct = aligned behavioral data by rat

high = cell(length(ratList), 1);
low = cell(length(ratList), 1);

for rat = 1:length(ratList)

    block = {'mixed'; 'high'; 'low'};
    T = linspace(-5, 10, size(pstruct.('CPIn'), 2)-1);
    bdata = bstruct.('CPIn');
    Rewards = convertreward(bdata.Reward);

    for bl = 2:length(block) %bl = block level

        %select trials
        if thirdrew_arg
            these = find(bdata.Block==bl & Rewards==3 ...
                & bdata.PrevTrialType~=2);
        else
            these = find(bdata.Block==bl ...
                & bdata.PrevTrialType~=2);
        end

        %get and save pstruct data with these criteria
        pdata = pstruct.('CPIn')(these, 2:end);
        da_mat = mean(pdata, 'omitnan');
        y = da_mat;
        if bl == 2
            bhi = y;
        else
            blo = y;
        end

    end

    %normalize by max for averaging across rats
    mymax = max([bhi blo], [], 'omitnan');
    high{rat, 1} = bhi/mymax;
    low{rat, 1} = blo/mymax;

end

end