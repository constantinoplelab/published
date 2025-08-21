function [high, low, high_err, low_err, T] = hilo_one_event(bstruct, ...
    pstruct, thirdrew_arg, event)
%%plot response aligned to event in high versus low blocks
%%input:
%ratlist = which rats to run through on at a time
%pstruct = aligned photometry data by rat
%bstruct = aligned behavioral data by rat

sem = @(xx) std(xx, 'omitnan') ./ sqrt(sum(~isnan(xx)));

block = {'mixed'; 'high'; 'low'};

T = linspace(-5, 10, size(pstruct.(event), 2)-1);

bdata = bstruct.(event);
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

    %get pstruct data with these criteria
    pdata = pstruct.(event)(these, 2:end);
    da_mat = mean(pdata, 'omitnan');

    %baseline correct
    y = da_mat;
    if bl == 2
        high = y;
        high_err = sem(pdata);
    else
        low = y;
        low_err = sem(pdata);
    end
end

end