function [xvec,spikes_binned] = binspikes(spiketimes, handles, alignment,dt,tmin,tmax,wndw,varargin)
%BINSPIKES put spikes into equalbins
%
%input:
%   spiketimes: long array of spike timings
%   
%   handles: handle of trial timings
%   
%   alignment: string handle for what event to be aligned
%   
%   dt: bin size (in s)
%   
%   tmin: first bin
%   
%   tmax: last bin
%   
%   wndw: time (in s) to extend from each side of tmin,tmax for spikes
%   if wndow if a length 2 vectyor, treat as asymmetric window for before
%   start and after end
%
%output:
%
%   xvec: tmin:dt:tmax. size nx
%   
%   spikes_pertrial: count of spikes. size (nt,nx)

if ~isempty(varargin)
    ops = varargin{1}; %helper structure to control binning
else
    ops = [];
end

if numel(wndw) ==2
    wndw_l = wndw(1);
    wndw_r = wndw(2);
else
    wndw_l = wndw;
    wndw_r = wndw;
end

ns = length(handles.start);
xvec = tmin:dt:tmax;
nt = numel(xvec);
spikes_binned = zeros(ns,nt);

%find trial spikes by searchign range between trial start and
%end/violation,
%+- wndw num seconds (for spikes, wndw= 1, stimuli wndw =0 are typical)

%find alignment data
switch alignment
    case 'start'
        alignHandle = handles.start;
    case 'choice'
        alignHandle = handles.choice;
    case 'leavecpoke'
        alignHandle = handles.leavecpoke;
    case 'end'
        alignHandle = handles.end;
    case 'lastflash'
        %go trial-wise and find
        %makes handling empty data eaasier than cellfun
        alignHandle = nan(size(handles.start));
        for j = 1:ns
            if isempty(handles.lflashes{j}) || isempty(handles.rflashes{j})
                continue
            else
                alignHandle(j) = max([handles.lflashes{j},handles.rflashes{j}]);
            end
        end
    otherwise %a custum alighment handle
        alignHandle = alignment;
end


%decide on handles.end or handles.leavec poke
%if rat participates in trial use handles.end
if isfield(handles,'end')
    endvec = handles.end;
else
    endvec = handles.leavecpoke;
end

%temporary boolean to turn on this non-overlap stuff
if isfield(ops,'restrictOverlap')
    restrictOverlap = ops.restrictOverlap;
else
    restrictOverlap = false;
end

for j = 1:ns

     %decide window size for each trial: winw_l/r, but prevent overlap with
    %previous and next trial
    if restrictOverlap
        if j > 1
            wndw_lj = min(wndw_l,handles.start(j)-endvec(j-1));
        else
            wndw_lj = wndw_l;
        end
        if j<ns
          wndw_rj = min(wndw_r,handles.start(j+1)-endvec(j));
        else
          wndw_rj = wndw_r;
        end
    else
        wndw_lj = wndw_l;
        wndw_rj = wndw_r;
    end

    these = spiketimes>=handles.start(j)-wndw_lj & spiketimes<=endvec(j)+wndw_rj;
    spikes_j = spiketimes(these)-alignHandle(j);
    
    
    %myu approach: sort into bins, aligned to a specific event   
    binedges = [xvec,xvec(end)+dt];
    yind = discretize(spikes_j,binedges);     
    %TODO: there is a non-loop version of doing this using ==. replace?
    for m = 1:nt
        spikes_binned(j,m) = sum(numel(spikes_j(yind==m)));
    end   
    
    %try using Christine's approach with hist
    %[n] = hist(spikes_j,xvec);
    %spikes_binned(j,:) = n;
    
end