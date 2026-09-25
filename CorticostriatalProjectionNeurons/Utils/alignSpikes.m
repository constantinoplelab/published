function [xvec, hmat, alignspikes] = alignSpikes(spiketimes, cond, wndw,...
    bins, alignto,trialend)
%This function takes spiketimes and aligns them to the event (cpoke,
%reward, opt out, etc) specified by alignto. Timestamps are according to 
%open ephys computer clock. cond is the condition found in
%behEvents (hits, high, low, etc). Generates a matrix called hmat where
%each row is a trial and each column is a timestep. The value is the firing
%rate in a bin specified by the input variable bins. 
% cmc 2019, SS 02/2020
%
%(desc by dlh, 09/2021)
%INPUTS:
%   spiketimes: array of spike times, in s
%   cond: (ns x 1) boolean of trial to keep. ns trials
%   wndw: (2 x1) array of times (in s) for keeping data. 
%       grabs [alignto + wndw(1), alignto + wndw(2)] window of data per trial
%   bins: time resolution (in s)
%   alignto: (nx x 1) array of timings of alignment event for each trial
%   trialend: (ns x 1) array of when each trial ends (in s).
%       each trial is trimmed to min([trialend, wndw(2) + alignto])
%
%OUTPUS:
%   xvec: (nt x 1) vector of time bins for firing firing rate
%   hmat: (ns x nt) heatmat of smoothed firing rates per trial
%   alignspikes: (ns x 1) cell of spike timings per trial, relative to
%       alignto
%

curr_trialend = trialend(:,2); %used as when rat pokes on following trial (so t+1)
prev_trialend = trialend(:,1); %used as when light comes on for current trial
%%time vector
xvec = [wndw(1):bins:wndw(2)];
d = diff(xvec)/2;
edges = [xvec(1)-d(1), xvec(1:end-1)+d, xvec(end)+d(end)];

if length(alignto) ~= length(cond)
    less = min([length(alignto) length(cond)]);
else
    less = length(alignto);
end

for j = 1:less
    t1 = alignto(j)+wndw(1);
    t2 = alignto(j)+wndw(2);

    %these are the spiketimes occurring in trial j during wndw
    these = find(spiketimes>=t1 & spiketimes<t2);
        
    %we want to make sure we are not including spikes from the next
    %trial.
    if j<length(alignto)
        nextT = curr_trialend(cond(j)+1);
        ia = find(spiketimes(these)>=nextT);
        these(ia) = [];
        beforeNextTrial = curr_trialend(cond(j)+1) - alignto(j);
    else
        nextT = [];
        beforeNextTrial = wndw(2);
    end
    %subtract align0
%     spiketimes(these) = spiketimes(these)-align0(cond(j));
    spiketimes(these) = spiketimes(these)-alignto(j);
    alignspikes{j,1} = spiketimes(these);
        
    [n,edges] = histcounts(spiketimes(these), edges);   
    
    hmat(j,:) = n./bins;
    hmat(j,:) = smooth(hmat(j,:));

    %NaN bins that bleed into next trial or come from previous trial
    if j > 1
        fromLastTrial = alignto(j) - prev_trialend(cond(j)); 
    else 
        fromLastTrial = abs(wndw(2));
    end
    
    if fromLastTrial < abs(wndw(1))
        time = abs(wndw(1)) - fromLastTrial;
        toNan = ceil(time/bins);

        if toNan > size(n,2)
             toNan = size(n,2);
        end

        xx = size(hmat,2);
        hmat(j,1:toNan) = nan;
        if xx~= size(hmat,2)
            keyboard
        end
        %hmat(:,length(xvec)+1:end)=[];  %added this
    end
    
    try
    if beforeNextTrial < wndw(2) 
        time = wndw(2)-beforeNextTrial;
        toNan = ceil(time/bins);
        if toNan > size(n,2)
             toNan = size(n,2);
        end
        hmat(j,length(xvec)-toNan:end) = nan;
    end
    catch
        keyboard
    end

    spiketimes(these) = spiketimes(these) + alignto(j);
   
end
 

end

