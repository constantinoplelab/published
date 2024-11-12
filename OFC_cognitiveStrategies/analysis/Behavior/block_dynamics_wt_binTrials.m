function [ltom, htom, mtol, mtoh, ltom_incong, htom_incong, sems] =...
    block_dynamics_wt_binTrials(A, twin, binSize, smoothfactor, varargin)
%computes the change in z-scored wait time, normalized to mean wait time 
%   for each reward volume around each block transition. **Does not include 
%   violation trials**. adapted from block_dynamics_wt. SSS 10/4/2023
%   inputs: A struct
%       twin = time window to analyze [-twin twin]
%       smoothfactor = smoothing factor. cmc has found 5-10 (trials) is 
%           reasonable for most rats.
%   outputs: a matrix of the average change in wait time for each
%            transition. 
%       ltom = low to mixed
%       htom = high to mixed
%       mtol = mixed to low
%       mtoh = mixed to high
  
if isfield(A, 'pd') %just in case the input is an S struct.
    [A, ~, ~] = parse_data_from_mysql(A);
end
[~, A.reward] = convertreward(A.reward);

A.wait_time(logical(A.catch)) = (A.wait_time(logical(A.catch))-...
    mean(A.wait_time(logical(A.catch)), 'omitnan'))./...
    std(A.wait_time(logical(A.catch)), 'omitnan');

%convert everything into z-scored wait time for that volume.
A = deltawt(A); 

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

for j = 1:length(lowtomix)
    these = (lowtomix(j)-twin):(lowtomix(j)+twin-1);
    [c, ia] = intersect(these, ctch);
    ltom(j,ia) =  A.wait_time(c);        
end

for j = 1:length(hightomix)
    these = (hightomix(j)-twin):(hightomix(j)+twin-1);
    [c, ia] = intersect(these, ctch);
    htom(j,ia) =  A.wait_time(c);
end

for j = 1:length(mixtolow)
    these = (mixtolow(j)-twin):(mixtolow(j)+twin-1);
    [c, ia] = intersect(these, ctch);
    mtol(j,ia) =  A.wait_time(c);
end

for j = 1:length(mixtohigh)
    these = (mixtohigh(j)-twin):(mixtohigh(j)+twin-1);
    [c, ia] = intersect(these, ctch);
    mtoh(j,ia) =  A.wait_time(c);
end

for j = 1:length(lowtomix_incong)
    these = (lowtomix_incong(j)-twin):(lowtomix_incong(j)+twin-1);
    [c, ia] = intersect(these, ctch);
    ltom_incong(j,ia) =  A.wait_time(c);
end

for j = 1:length(hightomix_incong)
    these = (hightomix_incong(j)-twin):(hightomix_incong(j)+twin-1);
    [c, ia] = intersect(these, ctch);
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


ltom = causal_smooth(ltom, smoothfactor);
htom = causal_smooth(htom, smoothfactor);
mtol = causal_smooth(mtol, smoothfactor);
mtoh = causal_smooth(mtoh, smoothfactor);

end

function [A] = deltawt(A)

rew = [5 10 20 40 80];
for j = 1:length(rew)
    m = find(A.reward==rew(j) & A.catch==1 & A.vios==0);
    A.wait_time(m) = (A.wait_time(m)-mean(A.wait_time(m), 'omitnan'))./...
        std(A.wait_time(m), 'omitnan');
end

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
    
end