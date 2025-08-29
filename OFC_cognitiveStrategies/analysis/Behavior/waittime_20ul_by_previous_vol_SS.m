function wait_time_mean = waittime_20ul_by_previous_vol_SS(A, inferredMixArg)
% Wait times for 20ul trials conditioned on whether the previous trial was
% <20 or >20

%   INPUTS: 
%       A - .mat file with rat data 
%   OUTPUT: 
%       wait_time_mean - Average wait time for < 20 uL and > 20 uL

[~, rew] = convertreward(A.reward);
wt = A.wait_time;

if ~isempty(inferredMixArg)
    BlkInf = inferredMixArg;
    find_20ul_catch = find(rew==20 & A.block==1 &...
        A.optout==1 & A.catch==1 & BlkInf==1); %20ul catch trials that are correctly inferred to be mixed block trials based on inference model
else
    find_20ul_catch = find(rew==20 & A.block==1 &...
        A.optout==1 & A.catch==1);
end

find_20ul_catch(find_20ul_catch == 1) = [];

wt = (wt - mean(wt(find_20ul_catch), 'omitnan'))./...
    std(wt(find_20ul_catch), 'omitnan');

less_than_20 = find_20ul_catch(rew(find_20ul_catch-1) < 20);
greater_than_20 = find_20ul_catch(rew(find_20ul_catch-1) > 20);

wait_time_mean = [mean(wt(less_than_20), 'omitnan'),...
    mean(wt(greater_than_20), 'omitnan')];

end

