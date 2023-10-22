function [ltom, htom, mtol, mtoh, ltom_incong, htom_incong, sems] =...
    block_dynamics_latency(A, twin, smoothfactor, varargin)
%plot how latency changes around block transitions. cmc 7/16/21
%inputs: A struct- If given S struct, will convert it to A.
%        twin- ntrials around block transitions.
%        smoothfactor- data points for smoothing.
%        varargin- whether to nan certain values
%outputs: ltom = mean delta z-scored latency at low to mixed transitions
%         htom = high to mixed
%         mtol = mixed to low
%         mtoh = mixed to high
%%there is also code below that, if uncommented, will generate synthetic
%%data with known kernels, for the purposes of debugging.

if isempty(varargin)
    detrendArg = true;
else
    detrendArg = varargin{1};
end

if isfield(A, 'pd') %just in case the input is an S struct.
    [A, ~, ~] = parse_data_from_mysql(A);
end
[~, A.reward] = convertreward(A.reward);

l = A.ITI;
l(cumsum([0; A.ntrials(1:end-1)])+1) = nan;
l(l>prctile(l, 99)) = nan;

ltncy = (l-mean(l, 'omitnan'))./std(l, 'omitnan');

if isempty(varargin)
    dropthese = true(size(ltncy));
else
    dropthese = varargin{1};
end
ltncy(~dropthese) = nan;

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

sems = nan(6, 2*twin+1);

[outmat] = meansub_latency(sess, treal, newy, twin, A, -2);

if ~isempty(outmat)
    ltom = mean(outmat, 1, 'omitnan');
    sems(1,:) = std(outmat, 'omitnan')./sqrt(size(outmat, 1));
    ltom(twin+1) = ltom(twin);
else
    ltom = nan(1, 2*twin+1);
    sems(1,:) = nan(1, 2*twin+1);
end

[outmat] = meansub_latency(sess, treal, newy, twin, A, -1);
if ~isempty(outmat)
    htom = mean(outmat, 1, 'omitnan');
    sems(2,:) = std(outmat, 'omitnan')./sqrt(size(outmat, 1));
    htom(twin+1) = htom(twin);
else
    htom = nan(1, 2*twin+1);
    sems(2,:) = nan(1, 2*twin+1);    
end

[outmat] = meansub_latency(sess, treal, newy, twin, A, 1);
if ~isempty(outmat)
    mtoh = mean(outmat, 1, 'omitnan');
    sems(3,:) = std(outmat, 'omitnan')./sqrt(size(outmat, 1));
    mtoh(twin+1) = mtoh(twin);
else
    mtoh = nan(1, 2*twin+1);
    sems(3,:) = nan(1, 2*twin+1);    
end

[outmat] = meansub_latency(sess, treal, newy, twin, A, 2);
if ~isempty(outmat)
    mtol = mean(outmat, 1, 'omitnan');
    sems(4,:) = std(outmat, 'omitnan')./sqrt(size(outmat, 1));
    mtol(twin+1) = mtol(twin-1);
else
    mtol = nan(1, 2*twin+1);
    sems(4,:) = nan(1, 2*twin+1);    
end

[outmat] = meansub_latency(sess, tincong, newy, twin, A, -2);
if ~isempty(outmat)
    ltom_incong = mean(outmat, 1, 'omitnan');
    sems(5,:) = std(outmat, 'omitnan')./sqrt(size(outmat, 1));
    ltom_incong(twin+1) = ltom_incong(twin);
else
    ltom_incong = nan(1, 2*twin+1);
    sems(5,:) = nan(1, 2*twin+1);    
end

[outmat] = meansub_latency(sess, tincong, newy, twin, A, -1);
if ~isempty(outmat)
    htom_incong = mean(outmat, 1, 'omitnan');
    sems(6,:) = std(outmat, 'omitnan')./sqrt(size(outmat, 1));
    htom_incong(twin+1) = htom_incong(twin);
else
    htom_incong = nan(1, 2*twin+1);
    sems(6,:) = nan(1, 2*twin+1);    
end

ltom(twin) = nan;
htom(twin) = nan;
mtol(twin) = nan;
mtoh(twin) = nan;

ltom = smooth(ltom, smoothfactor);
htom = smooth(htom, smoothfactor);
mtol = smooth(mtol, smoothfactor);
mtoh = smooth(mtoh, smoothfactor);

% xind = twin+1;
% ltom(xind:end) = smooth(ltom(xind:end), smoothfactor);
% ltom(1:twin) = flipud(smooth(fliplr(ltom(1:twin)), smoothfactor));
% 
% htom(xind:end) = smooth(htom(xind:end), smoothfactor);
% htom(1:twin) = flipud(smooth(fliplr(htom(1:twin)), smoothfactor));
% 
% mtol(xind:end) = smooth(mtol(xind:end), smoothfactor);
% mtol(1:twin) = flipud(smooth(fliplr(mtol(1:twin)), smoothfactor));
% 
% mtoh(xind:end) = smooth(mtoh(xind:end), smoothfactor);
% mtoh(1:twin) = flipud(smooth(fliplr(mtoh(1:twin)), smoothfactor));
% 
% ltom_incong(xind:end) = smooth(ltom_incong(xind:end), smoothfactor);
% ltom_incong(1:twin) = flipud(smooth(fliplr(ltom_incong(1:twin)),...
%     smoothfactor));
% 
% htom_incong(xind:end) = smooth(htom_incong(xind:end), smoothfactor);
% htom_incong(1:twin) = flipud(smooth(fliplr(htom_incong(1:twin)),...
%     smoothfactor));

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