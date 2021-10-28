function [width1mat,width2mat,wf_used,xplot,mask] = extractWaveforms(A,fileinds)
%extracts action potential (w1) and afterhyperpolarization (w2) widths

%% load fileinds and cluster ids

%G = load('/Users/dhocker/projects/ofc/data/published/concatdata_ofc.mat')

%neurons with troublesome waveforms. ones found in the first sweep
T = readtable('/Users/dhocker/projects/ofc/data/published/neurons_flagged.csv');
waveforms_flagged = T.neuron(T.keep==-1); %second pass

%% extract waveforms and look at spectral form

nn = numel(fileinds);

ns = numel(A{1}.waveform(:,1)); %number of entries
dt = 0.001; % sampling time, in s
T = ns*dt; %total time
Fmax = 1/dt; %max sampling freq, in Hz
Fmin = 1/T; %minimum sampling freq.

wvec = 0:Fmin:Fmax/2;
wave_w = zeros(nn,numel(wvec));

useinterp = true;

fs = 1/(32*1000);
xvec = 0:fs:fs*31;

if ~useinterp
    waveforms_all = zeros(nn,32,4); %not interpolated
else
    waveforms_all = zeros(nn,320*2,4); %interpolated
    xvecnew = linspace(0,xvec(end),320*2);
end

if useinterp
    xplot = xvecnew/1e-6;
else
    xplot = xvec/1e-6;
end

for j = 1:nn

    waveforms = A{fileinds(j)}.waveform;
    if ~useinterp
        waveforms_all(j,:,:) = waveforms;
    else
        waveforms_all(j,:,:) = interp1(xvec ,waveforms,xvecnew);
    end
    
    [~,idx_max] = max(max(waveforms,[],1));
    x = waveforms(:,idx_max)';

    %do fftr
    wave_wj = (1/(T/dt)) * abs(fft(x)).^2;
    wave_w(j,:) = wave_wj(1:numel(wvec));
   
end

%% extract width1 and width 2
mask = ~sum(fileinds == waveforms_flagged,1); %indices in fileinds that are usable
inds_used = find(mask); %index in fileinds for usable data. used for title plot
nn_used = numel(inds_used);

wf_used = waveforms_all(mask,:,:); %usable waveforms

%which channel to extract
[~,channel] = max(max(wf_used,[],2),[],3);
channel = squeeze(channel);

%location of max + peak
maxpeak = zeros(nn_used,1);
minpeak = zeros(nn_used,1);


for j = 1:nn_used
    [~,mp] = max(wf_used(j,:,channel(j)));
    maxpeak(j) = squeeze(mp);
    
    %min peak must occur after max peak
    [~,mp_inds] = sort(wf_used(j,:,channel(j)),'ascend');
    for k = 1:numel(mp_inds)
        if mp_inds(k) > maxpeak(j)
            minpeak(j) = mp_inds(k);
            break;
        end
    end
    
end

%find continuous regions above 0, continuous regions below 0
eps = 0.05;
width1mat = zeros(nn_used,2);
width2mat = zeros(nn_used,2);

for j = 1:nn_used
    %amplitude normalize 
    wf = wf_used(j,:,channel(j));
    r = wf(maxpeak(j))-wf(minpeak(j));
    waveform = (wf-mean(wf))/r;
    
    %find continous max region containing maxpeak
    posinds = find(waveform > eps);
    dpos = diff(posinds);
    breaks = [0,find(dpos ~= 1),numel(posinds)]; %location of region breaks

    
    %decide if that region is usable
    for k = 1:numel(breaks)-1
        indstart = posinds(breaks(k)+1);
        indstop = posinds(breaks(k+1));
        if indstart < maxpeak(j) && maxpeak(j) < indstop
            %width1 = [indstart, indstop];  %full peak width
            width1 = [maxpeak(j), indstop]; %half peak width
        end
    end
    
    %find continous max region containing maxpeak
    neginds = find(waveform < -eps);
    dneg = diff(neginds);
    breaks = [0,find(dneg ~= 1),numel(neginds)]; %location of region breaks
  
    %decide if that region is usable
    for k = 1:numel(breaks)-1
        
        indstart = neginds(breaks(k)+1);
        indstop = neginds(breaks(k+1));
        if indstart < minpeak(j) && minpeak(j) < indstop
            width2 = [indstart, indstop];           
        end
    end
    
    width1mat(j,:) = width1;
    width2mat(j,:) = width2;   
    
end









    
    