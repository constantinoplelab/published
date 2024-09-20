function [ltom, htom, mtol, mtoh,...
    ltom_incong, htom_incong, sems, outmat] =...
    block_dynamics_latency(A, twin, smoothfactor, z_arg, varargin)
% block_dynamics_latency - plot how latency changes around block 
% transitions
% INPUTS: 
%   A struct- Rat behavioral data. Usually BestFitEarly.(RAT_ID).All.ratTrial
%   twin- ntrials around block transitions.
%   smoothfactor- data points for smoothing.
%   varargin- three options:
%       1. Whether to detrend data or not (default true)
%       2. Whether there are certain trials you want to drop (default
%           none)
%       3. Whether to use causal smoothing (default false)
%           E.g., to change 3. leave 1. and 2. as default, use
%           block_dynamics_latency(A, twin, smoothfactor, z_arg,...
%               [], [], usecausal). That will use the default for 1. and 2.
% OUTPUTS: 
%   ltom = mean delta z-scored latency at low to mixed transitions
%   htom = high to mixed
%   mtol = mixed to low
%   mtoh = mixed to high

% Pull input parameters
if isempty(varargin) || isempty(varargin{1})
    detrendArg = true;
else
    detrendArg = varargin{1};
end

if length(varargin) < 2 || isempty(varargin{2})
    dropthese = true(size(A.ITI));
else
    dropthese = varargin{2};
end

%decide if doing causal smoothing, or standard moving window average
if length(varargin) < 3 || isempty(varargin{3})
    usecausal = false;
else
    usecausal = varargin{3};
end

[~, A.reward] = convertreward(A.reward); % Convert rewards to 1:5

% PUll ITIs
l = A.ITI;
l(cumsum([0; A.ntrials(1:end-1)])+1) = nan; % nan out first trial

% z-score itis
if z_arg
    l(l>prctile(l, 99)) = nan;
    ltncy = (l-mean(l, 'omitnan'))./std(l, 'omitnan');
else
    ltncy = l;
end

% Nan out drop conditions
ltncy(~dropthese) = nan;

% Find block transitions
[~, treal, tincong] = blockhmat_withvios(A);

%initialize matrix for latencies
sess = nan(length(A.ntrials), max(A.ntrials)); 
ctr = 1;

xvec = (1:max(A.ntrials))';
s = []; x = [];
for jk = 1:length(A.ntrials)
    sess(jk,1:A.ntrials(jk)) = ltncy(ctr:ctr+A.ntrials(jk)-1);
    x = [x; xvec];
    s = [s; sess(jk,:)'];
    ctr = ctr+A.ntrials(jk);
end

%%we want to regress latency against trial number, so we can subtract
%%it out later.
if detrendArg
    bad = find(isnan(s));
    s(bad) = [];
    x(bad) = [];
    X = [x, ones(length(x),1)];
    [b] = regress(s,X);
    newy = (b(1).*xvec + b(2))'; %average latency by trial number.
else
    newy = zeros(1, size(sess, 2));
end

% Preallocate data structures
sems = nan(4, 2*twin+1);
outmat = cell(4, 1);

% Pull low to mixed
[outmat{1}] = meansub_latency(sess, treal, newy, twin, A, -2);
if ~isempty(outmat{1})
    ltom = mean(outmat{1}, 1, 'omitnan');
    sems(1,:) = std(outmat{1}, 'omitnan')./sqrt(size(outmat{1}, 1));
    ltom(twin+1) = ltom(twin);
else
    ltom = nan(1, 2*twin+1);
    sems(1,:) = nan(1, 2*twin+1);
end

% Pull high to mixed
[outmat{2}] = meansub_latency(sess, treal, newy, twin, A, -1);
if ~isempty(outmat{2})
    htom = mean(outmat{2}, 1, 'omitnan');
    sems(2,:) = std(outmat{2}, 'omitnan')./sqrt(size(outmat{2}, 1));
    htom(twin+1) = htom(twin);
else
    htom = nan(1, 2*twin+1);
    sems(2,:) = nan(1, 2*twin+1);    
end

% Pull mixed to high
[outmat{3}] = meansub_latency(sess, treal, newy, twin, A, 1);
if ~isempty(outmat{3})
    mtoh = mean(outmat{3}, 1, 'omitnan');
    sems(3,:) = std(outmat{3}, 'omitnan')./sqrt(size(outmat{3}, 1));
    mtoh(twin+1) = mtoh(twin);
else
    mtoh = nan(1, 2*twin+1);
    sems(3,:) = nan(1, 2*twin+1);    
end

% Pull mixed to low
[outmat{4}] = meansub_latency(sess, treal, newy, twin, A, 2);
if ~isempty(outmat{4})
    mtol = mean(outmat{4}, 1, 'omitnan');
    sems(4,:) = std(outmat{4}, 'omitnan')./sqrt(size(outmat{4}, 1));
    mtol(twin+1) = mtol(twin-1);
else
    mtol = nan(1, 2*twin+1);
    sems(4,:) = nan(1, 2*twin+1);    
end

% Smooth data
ltom(twin) = nan;
htom(twin) = nan;
mtol(twin) = nan;
mtoh(twin) = nan;

if usecausal
    ltom = causal_smooth(ltom, smoothfactor);
    htom = causal_smooth(htom, smoothfactor);
    mtol = causal_smooth(mtol, smoothfactor);
    mtoh = causal_smooth(mtoh, smoothfactor);
else
    ltom = smooth(ltom, smoothfactor);
    htom = smooth(htom, smoothfactor);
    mtol = smooth(mtol, smoothfactor);
    mtoh = smooth(mtoh, smoothfactor);
end

end


function [outmat] =...
    meansub_latency(latencymat, diffmat, regline, wndw, A, trialtype)

ctr = 1;
outmat = nan(1, 2*wndw+1);
for xm = 1:length(diffmat(:,1))
    th = find(diffmat(xm,:)==trialtype);
    for jk = 1:length(th)
        these = th(jk)-wndw:th(jk)+wndw;        
        these(these>length(latencymat(1,:))) = [];
        
        [tt, ia] = intersect(these, 1:length(A.reward));
        
        outmat(ctr,ia) = latencymat(xm,tt) - regline(1,tt);
        ctr = ctr+1;
    end
end

end