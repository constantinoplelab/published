function [ltom, htom, mtol, mtoh, ltom_incong, htom_incong, sems] =...
    block_dynamics_wt_noNormalization(A, twin, binSize, smoothfactor, ...
    volType, varargin)
%computes the change in raw wait time around each block transition. **Does not include 
%   violation trials**. adapted from block_dynamics_wt_binTrials. SSS 2025

%   INPUTS: 
%       A struct: .mat file with behavior data for individual rat
%       twin = time window to analyze [-twin twin]
%       binSize = how many bins to average over
%       smoothfactor = smoothing factor. cmc has found 5-10 (trials) is 
%           reasonable for most rats.
%       volType = which volumes to include in analysis (1 = congruent
%           volumes; 0 = 20ul only)
%       varargin = indices to look at specific trial sets

%   OUTPUTS: a matrix of the average change in wait time for each
%            transition. 
%       ltom = low to mixed
%       htom = high to mixed
%       mtol = mixed to low
%       mtoh = mixed to high
  

[~, A.reward] = convertreward(A.reward);

if ~isempty(varargin)
    [~, ia] = setdiff(1:length(A.wait_time), varargin{1}); %set elements not in varargin to nan
    A.wait_time(ia) = nan;
end


A = removeVios(A); %remove violation trials


[lowtomix, hightomix,...
    mixtolow, mixtohigh,...
    lowtomix_incong, hightomix_incong] =...
    find_block_transition(A);

xvec = -twin:5:twin-1;
ltom = nan(length(lowtomix), twin*2);
htom = nan(length(hightomix), twin*2);
mtol = nan(length(mixtolow), twin*2);
mtoh = nan(length(mixtohigh), twin*2);
ltom_incong = nan(length(lowtomix), twin*2);
htom_incong = nan(length(hightomix), twin*2);

ctch = find(A.catch==1);

if volType == 1
    % use all congruent volume trials for each transition type
    large_catch = intersect(intersect(find(A.block == 1 | A.block == 2), ...
        find(A.reward == 40 | A.reward == 80)), ctch);
    small_catch = intersect(intersect(find(A.block == 1 | A.block == 3), ...
        find(A.reward == 5 | A.reward == 10)), ctch);

elseif volType == 0
    % only use 20ul trials
    large_catch = intersect(intersect(find(A.block == 1 | A.block == 2), ...
        find(A.reward == 20)), find(A.catch==1));
    small_catch = intersect(intersect(find(A.block == 1 | A.block == 3), ...
        find(A.reward == 20)), find(A.catch==1));
end


for j = 1:length(lowtomix)
    these = (lowtomix(j)-twin):(lowtomix(j)+twin-1);
    [c, ia] = intersect(these, small_catch);
    ltom(j,ia) =  A.wait_time(c);        
end

for j = 1:length(hightomix)
    these = (hightomix(j)-twin):(hightomix(j)+twin-1);
    [c, ia] = intersect(these, large_catch);
    htom(j,ia) =  A.wait_time(c);
end

for j = 1:length(mixtolow)
    these = (mixtolow(j)-twin):(mixtolow(j)+twin-1);
    [c, ia] = intersect(these, small_catch);
    mtol(j,ia) =  A.wait_time(c);
end

for j = 1:length(mixtohigh)
    these = (mixtohigh(j)-twin):(mixtohigh(j)+twin-1);
    [c, ia] = intersect(these, large_catch);
    mtoh(j,ia) =  A.wait_time(c);
end

for j = 1:length(lowtomix_incong)
    these = (lowtomix_incong(j)-twin):(lowtomix_incong(j)+twin-1);
    [c, ia] = intersect(these, small_catch);
    ltom_incong(j,ia) =  A.wait_time(c);
end

for j = 1:length(hightomix_incong)
    these = (hightomix_incong(j)-twin):(hightomix_incong(j)+twin-1);
    [c, ia] = intersect(these, large_catch);
    htom_incong(j,ia) =  A.wait_time(c);
end

sems = nan(6, 2*twin/binSize);

sem = @(x) std(x, [], 'omitnan')./ sqrt(sum(~isnan(x))); %exclude nans

ltom = reshape(ltom, [], twin*2/binSize);
sems(1,:) = sem(ltom);
ltom = mean(ltom, 1, 'omitnan'); 

htom = reshape(htom, [], twin*2/binSize);
sems(2,:) = sem(htom);
htom = mean(htom, 1, 'omitnan');

mtol = reshape(mtol, [], twin*2/binSize);
sems(3,:) = sem(mtol);
mtol = mean(mtol, 1, 'omitnan');

mtoh = reshape(mtoh, [], twin*2/binSize);
sems(4,:) = sem(mtoh);
mtoh = mean(mtoh, 1, 'omitnan');

ltom_incong = reshape(ltom_incong, [], twin*2/binSize);
sems(5,:) = sem(ltom_incong);
ltom_incong = mean(ltom_incong, 1, 'omitnan');

htom_incong = reshape(htom_incong, [], twin*2/binSize);
sems(6,:) = sem(htom_incong);
htom_incong = mean(htom_incong, 1, 'omitnan');


ltom = causal_smooth_SS(ltom, smoothfactor);
htom = causal_smooth_SS(htom, smoothfactor);
mtol = causal_smooth_SS(mtol, smoothfactor);
mtoh = causal_smooth_SS(mtoh, smoothfactor);

end

function [lowtomix, hightomix,...
    mixtolow, mixtohigh,...
    lowtomix_incong, hightomix_incong] =...
    find_block_transition(A)

ntrials = [0; cumsum(A.ntrials)];

[~, treal, tincong] = blockhmat_withvios(A);

lowtomix = cell(length(ntrials), 1);
hightomix = cell(length(ntrials), 1);
mixtolow = cell(length(ntrials), 1);
mixtohigh = cell(length(ntrials), 1);
lowtomix_incong = cell(length(ntrials), 1);
hightomix_incong = cell(length(ntrials), 1);

for sess = 1:length(ntrials)-1
    lowtomix{sess} = (find(treal(sess,:) == -2) + ntrials(sess))';
    hightomix{sess} = (find(treal(sess,:) == -1) + ntrials(sess))';
    
    mixtolow{sess} = (find(treal(sess,:) == 2) + ntrials(sess))';
    mixtohigh{sess} = (find(treal(sess,:) == 1) + ntrials(sess))';
    
    lowtomix_incong{sess} = (find(tincong(sess,:) == -2) +...
        ntrials(sess))';
    hightomix_incong{sess} = (find(tincong(sess,:) == -1) +...
        ntrials(sess))';
end

lowtomix = cell2mat(lowtomix);
hightomix = cell2mat(hightomix);
mixtolow = cell2mat(mixtolow);
mixtohigh = cell2mat(mixtohigh);
lowtomix_incong = cell2mat(lowtomix_incong);
hightomix_incong = cell2mat(hightomix_incong);

end

function [A] = removeVios(A)
    vec = [0; A.ntrials];

    ntrials = arrayfun(@(x) A.ntrials(x) - ...
        sum(A.vios(sum(vec(1:x))+1:sum(vec(1:x+1)))), 1:length(vec)-1)';
   
    A.ntrials = ntrials;

    A.wait_time = A.wait_time(A.vios == 0);
    A.reward = A.reward(A.vios == 0);
    A.block = A.block(A.vios == 0);
    A.catch = A.catch(A.vios==0);
    A.optout = A.optout(A.vios == 0);
    
    
end