function [pvalue, permute_dist] = permute_test(a, b, numsamps, varargin)
%non-parametric permutation test for comparing values in a and b (a and b are column vectors).
%numsamps is the number of permutations. 100 or 1000 are reasonable. cmc 5/13/20

% a should still be the data with fewer samples
if isempty(varargin)
    data = abs(mean(a, 'omitnan')-mean(b, 'omitnan')); %use absolute value so order of a and b doesn't matter
else
    data = mean(a, 'omitnan') - mean(b, 'omitnan');
end
mix = [a; b];
new_a = [];
new_b = [];


for j = 1:numsamps
  samp_indx = ceil(rand(length(a), 1) * length(mix));
  new_a = [new_a; mean(mix(samp_indx), 'omitnan')];
  samp_indx = setdiff(1:length(mix), samp_indx);
  new_b = [new_b; mean(mix(samp_indx), 'omitnan')];
end

%null distribution
permute_dist = new_a-new_b; 

%integrate null dist. at data value
pvalue = numel(find(permute_dist>data))./numel(permute_dist); 
