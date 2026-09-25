function [hi, lo, mix] = blocks(A, varargin)
%script to calculate wait times in each block type. 
%inputs: A struct or S struct (which will then be converted to A struct).
%        varargin, if it exists, specifies trial indices to analyze
%output: hi, lo, and mix structs with fields 'wt' and 'er', corresponding to mean
%wait time and s.e.m. for each volume in each block.
%modified by cmc 10/19/21

%block 1 is test block.
%block 2 = high adaptation block.
%block 3 = low adaptation block.

if isfield(A, 'pd') %just in case the input is an S struct.
    [A, ~, ~] = parse_data_from_mysql(A);
end

rew = 1:5;
A.reward = convertreward(A.reward);

hi.wt = nan(1, length(rew));
hi.er = nan(1, length(rew));

lo.wt = nan(1, length(rew));
lo.er = nan(1, length(rew));

mix.wt = nan(1, length(rew));
mix.er = nan(1, length(rew));


for j = 1:length(rew)
    
    ix = find(A.reward==rew(j) & A.block==1 & A.optout==1);
    if ~isempty(varargin)
        ix = intersect(ix, varargin{1});
    end
    mix.wt(1,j) = mean(A.wait_time(ix), 'omitnan');
    mix.er(1,j) = std(A.wait_time(ix), 'omitnan')./sqrt(length(ix));
    
    ix = find(A.reward==rew(j) & A.block==2 & A.optout==1);
    if ~isempty(varargin)
        ix = intersect(ix, varargin{1});
    end
    hi.wt(1,j) = mean(A.wait_time(ix), 'omitnan');
    hi.er(1,j) = std(A.wait_time(ix), 'omitnan')./sqrt(length(ix));
    
    ix = find(A.reward==rew(j) & A.block==3 & A.optout==1);
    if ~isempty(varargin)
        ix = intersect(ix, varargin{1});
    end
    lo.wt(1,j) = mean(A.wait_time(ix), 'omitnan');
    lo.er(1,j) = std(A.wait_time(ix), 'omitnan')./sqrt(length(ix));
end

L  = A.ITI;
L(L>prctile(L,99)) = NaN;
hi.lat = mean(L(A.block==2), 'omitnan');
hi.lat_er = std(L(A.block==2), 'omitnan')./...
    sqrt(sum(A.block==2 & ~isnan(L)));

lo.lat = mean(L(A.block==3), 'omitnan');
lo.lat_er = std(L(A.block==3), 'omitnan')./...
    sqrt(sum(A.block==3 & ~isnan(L)));