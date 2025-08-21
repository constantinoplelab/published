function [DA,  T, tzero, delay_legend, bins] =...
    delay_effect_cg(bstruct, pstruct, event,...
    block, delaybins)
%%Plot dF/F as a function of reward delay. Exclude first 2s following side
%%LED to avoid cue response (see SideLED_resp.m for reference) cmc 8/6/20.
%currently, it's using screened behavioral data 9/1/21 CG
%t3bl = -1
%site = which site for the naming of the file (1 or 2)

% define bins
bins = delaybins;

%initiate variables
DA = cell(1, length(bins));
T = linspace(-5, 10, size(pstruct.(event), 2)-1);
[~ , tzero] = min(abs(T+0));
bdata = bstruct.(event);

%index trial types
these = find(bdata.Block==block & bdata.Catch==0 &...
    bdata.OptOut==0);        
pdata = pstruct.(event)(these, 2:end);

%find delay
delay = bdata.RewardDelay(these);
delay_legend = cell(length(bins), 1);
for ka = 1:length(bins)
    delay_legend{ka} = num2str(bins(ka));
end
d = delay;

delta = nan(1,length(bins));
for j = 2:length(delta)
    ix = d>bins(j-1) & d<=bins(j);
    da_mat = mean(pdata(ix,:), 'omitnan');
    DA{1, j} = da_mat;
end

end
