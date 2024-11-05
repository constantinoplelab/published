function [times_resample, sorted_data] = heatmat_byblock(pstruct, bstruct)

resample_freq = 0.1;
T = linspace(-5, 10, size(pstruct.('CPIn'), 2)-1);
bdata = bstruct.('CPIn');
rewards = convertreward(bdata.Reward);
blocks = bdata.Block;

%reorder blocks into low=1, mixed=2, high=3
blocks(blocks==1) = 0;
blocks(blocks==3) = 1;
blocks(blocks==2) = 3;
blocks(blocks==0) = 2;

%select trials
these = rewards==3 & blocks~=2 & bdata.PrevTrialType~=2; %only use 16 ul, adaptation blocks, and non-post-violation trials
pdata = pstruct.('CPIn')(these, 2:end);
blocks_these = blocks(these);

%down-sample
times_resample = -5:resample_freq:10; 
resample_indx = discretize(T, times_resample);
da_mat = arrayfun(@(ii) mean(pdata(:, resample_indx == ii), 2,...
    'omitnan'), unique(resample_indx), 'UniformOutput', 0);
da_mat = cell2mat(da_mat); 
[~ , tzero] = min(abs(times_resample+0));
[~ , tafter] = min(abs(times_resample-0.5));

%sort
maxpeak = max(da_mat(:,tzero:tafter),[],2);
A = [da_mat, blocks_these, maxpeak];
sorted_data = sortrows(A,[size(A,2)-1 size(A,2)]);

end

