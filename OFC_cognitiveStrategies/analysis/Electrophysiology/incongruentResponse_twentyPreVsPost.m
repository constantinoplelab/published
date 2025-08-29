function [pre20_pref, pre20_nonpref, post20_pref, post20_nonpref, sig] = ...
    incongruentResponse_twentyPreVsPost(SU, S, alignto, zarg)
% Compute responses on the last 20ul trial before and the first 20ul trial
% after an incongruent trial. Named based on transition type (i.e. a
% preferred transition is low to mix for a block that fires at a higher
% rate in high blocks)


%get incongruent trials
ntrials = length(S.Block);
io = getIncongruentTrials(S); %find the indices of the first incongruent offer after block transition
bad = find(io.ltom == ntrials); %no post-incongruent trials if last trial
io.ltom(bad,:) = [];
bad = find(io.htom == ntrials);
io.htom(bad,:) = [];


%convert volumes
v = convertreward(S.RewardAmount);

%responses in high and low blocks on average
n = length(SU);
w = arrayfun(@(x) find(round(SU(1).xvec.COFF, 3) == x), [0.25 0.5]); %using [0 0.5] as wndw - to calculate block significance

zSU = arrayfun(@(x) SU(x).hmat.(alignto), 1:n, 'uniformoutput', false); %raw fr

byBlock = arrayfun(@(y) arrayfun(@(x) mean(zSU{y}(S.Block == x, ...
    w(1):w(2)), 2, 'omitnan'), [2, 3], 'UniformOutput', false), 1:n, ...
    'uniformoutput', false); %base on raw firing rates

means = arrayfun(@(y) arrayfun(@(x) mean(byBlock{y}{x}, 'omitnan'), 1:2), 1:n, ...
    'uniformoutput', false); %means(1) = high, means(2) = low

%test for significant difference
h = arrayfun(@(x) ttest2(byBlock{x}{1}, byBlock{x}{2}), 1:n);

%specifically look at cells with significant block/volume encoding
zSU = zSU(h==1); 
means = means(h==1);
nc = length(zSU);

%keep track of which cells are significant
sig = zeros(length(SU), 1);
sig(h==1) = 1;

%get preferences for each cell
vec = 2./ones(nc,1);
prefLo = find(arrayfun(@(x) means{x}(1) < means{x}(2), 1:nc)); %cells that fire at a higher rate in low blocks
prefvec = vec; prefvec(prefLo) = 3;


%mixed 20 non-violation trials
block20 = arrayfun(@(x) intersect(intersect(find(S.Block == x), ...
    find(S.vios == 0)), find(v == 3)), 1:3, 'UniformOutput', false);


%high to mix transitions
nt = size(io.htom, 1);
post20_htom = cell2mat(arrayfun(@(x) block20{1}(intersect(find(block20{1} > ...
    io.htom(x, 1), 1), find(block20{1} < (io.htom(x)+10)))), 1:nt, ...
    'UniformOutput', false)'); %first rewarded 20ul trial post incongruent trial


himix = [block20{1}; block20{2}];
pre20_htom = cell2mat(arrayfun(@(x) himix(intersect(find(himix < io.htom(x,1), 1, 'last'), ...
    find(himix > (io.htom(x,1)-10)))), 1:nt, 'UniformOutput', false)'); %only the last trial pre incongruent

%low to mix transitions
nt = size(io.ltom, 1);
post20_ltom = cell2mat(arrayfun(@(x) block20{1}(intersect(find(block20{1} > ...
    io.ltom(x, 1), 1), find(block20{1} < (io.ltom(x)+10)))), 1:nt, ...
    'UniformOutput', false)');

lomix = [block20{1}; block20{3}];
pre20_ltom = cell2mat(arrayfun(@(x) lomix(intersect(find(lomix < io.ltom(x,1), 1, 'last'), ...
    find(lomix > (io.ltom(x,1)-10)))), 1:nt, 'UniformOutput', false)'); %only the last trial pre incongruent


%responses
if zarg
    zscore = @(x) (x - mean(x(:),'omitnan')) / std(x,0,'all','omitnan');

    zSU = arrayfun(@(x) zscore(zSU{x}), 1:nc, 'uniformoutput', ...
        false);
end

check_ltom = ~isempty(pre20_ltom) && ~isempty(post20_ltom);
check_htom = ~isempty(pre20_htom) && ~isempty(post20_htom);

%pre-allocate matrices based on all cells. Only cells with significant
%   difference in firing rate will be filled in
pre20_pref = nan(n, size(SU(1).hmat.COFF, 2));
pre20_nonpref = pre20_pref;
post20_pref = pre20_pref;
post20_nonpref = pre20_pref;

for ii = 1:nc
    a = find(h);
    idx = a(ii); %add into correct index based on all cells
    if prefvec(ii) == 2 
        if check_htom
            pre20_nonpref(idx,:) = mean(zSU{ii}(pre20_htom, :), 1, 'omitnan'); %non-preferred transition
            post20_nonpref(idx,:) = mean(zSU{ii}(post20_htom, :), 1, 'omitnan');
        end

        if check_ltom
            pre20_pref(idx,:) = mean(zSU{ii}(pre20_ltom, :), 1, 'omitnan'); %preferred transition
            post20_pref(idx,:) = mean(zSU{ii}(post20_ltom, :), 1, 'omitnan');
        end

    else
        if check_htom
            pre20_pref(idx,:) = mean(zSU{ii}(pre20_htom, :), 1, 'omitnan'); %preferred transition 
            post20_pref(idx,:) = mean(zSU{ii}(post20_htom, :), 1, 'omitnan');
        end

        if check_ltom
            pre20_nonpref(idx,:) = mean(zSU{ii}(pre20_ltom, :), 1, 'omitnan'); %non-preferred transition
            post20_nonpref(idx,:) = mean(zSU{ii}(post20_ltom, :), 1, 'omitnan');
        end
    end

end

