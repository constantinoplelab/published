function [times_resample, sorted_data] = heatmat_byreward(pstruct, bstruct)

resample_freq = 0.1;
            
T = linspace(-5, 10, size(pstruct.('CPIn'), 2)-1);

bdata = bstruct.('CPIn');
rewards = convertreward(bdata.Reward);

%select trials
these = bdata.Block==1 & bdata.PrevTrialType~=2; 
pdata = pstruct.('CPIn')(these, 2:end);
rewards_these = rewards(these);

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
A = [da_mat, rewards_these, maxpeak];
sorted_data = sortrows(A,[size(A,2)-1 size(A,2)]);

end

