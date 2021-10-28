%script to compress the start, choice, and poke data so that you aren't opening up tons of files

cd('~/projects/ofc/data/PhysiologyData/Christine_OFC_2019')
%savedir = '~/projects/ofc/data/PhysiologyData/Christine_OFC_2019';
[fnames, units, depth, duplicate, rat, days] = getfnames();


savetype = 'published';
switch savetype
    case 'published'
        savedir = '~/projects/ofc/data/published/'; 
        %only include usable neurons: > 1Hz single units and had decent neuron fits. 
        fileinds_name = '/Users/dhocker/projects/ofc/results/fileinds_659.mat';
        f = load(fileinds_name);
        fileinds = f.fileinds;
        fnames = fnames(fileinds);
    case 'processed'
        savedir = '~/projects/ofc/data/processed/'; %for my data
end

%if saving to different location


fname_save = strcat(savedir,'concatdata_ofc.mat'); %saved file

longForm = false; %boolean for whether or not saving long format, which includes clicks and flashes


%current structure
%A: cell array indexed by file order from fnames
%A{i} is a struc object with following fields
%   fname: the file name of the raw data
%   hmat*: the binned firing rates, aligned to that event
%   xvec*: the timing of the bins in hmat, aligned ot that event
%   handles: the explicit timing information for each trial
%   hits: if each trial was a success or not (1=yes,0=no,nan=violated
%   trial)
%   prev_hits: if previous trial was success or not
%   (1=yes,0=no,nan=violated)
%   went_right: was the right-side port chosen?
%   left_vol/right_vol: left port/right port volume of water on each trial
%   left_prob/right_prob: left and right port probabilities
%   chosenprob: which probability was chosen on that trial
%   chosenval: which volume of water was chosen on that trial
%   nspikes: total number of spikes per trial
%   fname: files name data taken from
%   isUsable: boolean for whether or not to include in analysis. currently
%       exluces multi-unit, and very low firing rate
%   selectivity: results of t-tests for whether or not number of spikes was
%       significantly higher on trials in a particular condition. +1 = case 1
%       had more spikes, -1 = case 2 had more spikes, 0 = nonselective
%   ratname: which animal was used.
%   trial_idx: trial number
%   depth: estimated depth (mm) of tetrode during recording
%   waveform: average waveforms on the 4 tetrode channels
%   unitType: single or multiunit
%   duplicate: was another recorded neuron form a separate session


%create a data structure
A = cell(numel(fnames),1);

for j = 1:numel(fnames)
    if mod(j,100) == 0
        disp(j)
    end
    
    load(strcat(['parsed_data', filesep, fnames{j}, '.mat'])); %hmats were generated using save_heatmats.m
   
    %run the align to leave cpoke
    [hmat_leavecpoke,xvec_leavecpoke] = alignhmattoleavecpoke(S,spiketimes,4,0.05,handles);

    [chosenprob, chosenval, hits, ~] = parse_choices(S);
    prev_hits = [nan;hits(1:end-1)];
    
    switch savetype
        case 'published'
            A{j}.fname = strcat('unit_',num2str(j),'.mat');
        case 'processed'
            A{j}.fname = fnames{j};
    end
    A{j}.fileind = j;

    %add the heat mats aligned to different events. save xvec outside to
    %avoid redundancy

    %rebin the hmats, rather than using ones from files. removes
    %odd edge behavior, but BE SURE to tell Chritine about this
    dt = 0.05;
   
    %start
    [~,hmat_start_rebin] = binspikes(spiketimes, handles, 'start',dt,-2,4,[2,4]);
    for k = 1:size(hmat_start_rebin,1)
        hmat_start_rebin(k,:) = smooth(hmat_start_rebin(k,:))/dt;
    end
    A{j}.hmat_start = hmat_start_rebin;
    %A{j}.hmat_start = hmat_start; old form
   
    %pokeend. don't know what this correspond to...
    A{j}.hmat_pokeend = hmat_pokeend; 
   
    %choice
    [~,hmat_choice_rebin] = binspikes(spiketimes, handles, 'choice',dt,-4,4,[4,4]);
    for k = 1:size(hmat_choice_rebin,1)
    hmat_choice_rebin(k,:) = smooth(hmat_choice_rebin(k,:))/dt;
    end
    A{j}.hmat_choice = hmat_choice_rebin;
    %A{j}.hmat_choice = hmat_choice; %old ofmr
   
    %lastflash
    [~,hmat_lastflash_rebin] = binspikes(spiketimes, handles, 'lastflash',dt,-4,4,[4,4]);
    for k = 1:size(hmat_lastflash_rebin,1)
        hmat_lastflash_rebin(k,:) = smooth(hmat_lastflash_rebin(k,:))/dt;
    end
    A{j}.hmat_lastflash = hmat_lastflash_rebin;
    %A{j}.hmat_lastflash = hmat_lastflash;
   
    %end
    [~,hmat_end_rebin] = binspikes(spiketimes, handles, 'end',dt,-4,4,[4,4]);
    for k = 1:size(hmat_end_rebin,1)
        hmat_end_rebin(k,:) = smooth(hmat_end_rebin(k,:))/dt;
    end
    A{j}.hmat_end = hmat_end_rebin;
    %A{j}.hmat_end = hmat_end;
   
    %leavecpoke (my favorite)
    [~,hmat_leavecpoke_rebin] = binspikes(spiketimes, handles, 'leavecpoke',dt,-4,4,[4,4]);
    for k = 1:size(hmat_leavecpoke_rebin,1)
        hmat_leavecpoke_rebin(k,:) = smooth(hmat_leavecpoke_rebin(k,:))/dt;
    end
    A{j}.hmat_leavecpoke = hmat_leavecpoke_rebin;
    %A{j}.hmat_leavecpoke = hmat_leavecpoke;
   
    %add timing information, with flashes and beeps moved outside of struct
    handles_binned = struct();
    handles_binned.start = handles.start;
    handles_binned.end = handles.end;
    %handles_binned.went_right = handles.went_right;
    handles_binned.choice = handles.choice;
    handles_binned.leavecpoke = handles.leavecpoke;
   
    %should lflashes and bups timing be included???
    if longForm
        handles_binned.lflashes = handles.lflashes;
        handles_binned.rflashes = handles.rflashes;
        handles_binned.lbups = handles.lbups;
        handles_binned.rbups = handles.rbups;
    end
   
    
    A{j}.handles = handles_binned;
    
    A{j}.hits = hits;
    A{j}.prev_hits = prev_hits;
    A{j}.went_right = handles.went_right; %different than parse_choice
    A{j}.this_left_volume = S.pd{:}.this_left_volume;
    A{j}.this_right_volume = S.pd{:}.this_right_volume;
    A{j}.left_prob = S.pd{:}.left_prob;
    A{j}.right_prob = S.pd{:}.right_prob;
    A{j}.chosenprob = chosenprob;
    A{j}.chosenval = chosenval;
   
    %classify cells as useable in analysis,
    %having  >=2 spikes on half of trials
    n = nspikespertrials(spiketimes, handles, 1);
    A{j}.nspikes = n;
    
    
    A{j}.ratname = S.ratname;
   
   
    A{j}.trial_idx = 1:size(A{j}.hmat_start,1);
    A{j}.waveform = waveform;
    
    
    %add selectivity information for processed data only. omit 4 publisehd
    %+1 means higher rate for first case, -1 higher for second case. 0 nonselective 
    %also include duplicate and multi vs single
    switch savetype
        case 'processed'
            
            %is usable beause it had enough spikes?
            nk = n>=2;
            if nanmean(nk) < 0.5 
                A{j}.isUsable = false;
            elseif ~strcmp(units{j},'single') %check if single unit or multi
                A{j}.isUsable = false;
            else
                A{j}.isUsable = true;
            end
            
            %selectivity?
            [n_LR,n_winloss, n_prewinloss, n_safeRisky,~,~,~,~] = ofc_Ttests(A(j));
            selectivity = struct();
            selectivity.LeftRight = n_LR;
            selectivity.WinLoss = n_winloss;
            selectivity.preWinLoss = n_prewinloss;
            selectivity.SafeRisky = n_safeRisky;

            A{j}.selectivity = selectivity; 
            
            A{j}.duplicate = duplicate;
            
            %there are some typos in units. check
            if contains(units(j),'single')
                A{j}.unitType = 'single';
            elseif contains(units(j),'multi')
                A{j}.unitType = 'single';
            else
                 A{j}.unitType = 'single';
            end
            
            A{j}.depth = depth(j);
    end
    
    
    
    

end

disp('saving')
save(fname_save,'A','xvec_*','-v7.3')

%%   
%add binned version of flahses and bups
   %{
   [ns,nt_start] = size(hmat_start);
   [~,nt_lastflash] = size(hmat_lastflash);
   [~,nt_choice] = size(hmat_choice);
   [~,nt_end] = size(hmat_end);
   [~,nt_pokeend] = size(hmat_pokeend);
   
   lflashCt_start = zeros(ns,nt_start);
   lflashCt_lastflash = zeros(ns,nt_lastflash);
   lflashCt_choice = zeros(ns,nt_choice);
   lflashCt_end = zeros(ns,nt_end);
   lflashCt_pokeend = zeros(ns,nt_pokeend);
   
   rflashCt_start = zeros(ns,nt_start);
   rflashCt_lastflash = zeros(ns,nt_lastflash);
   rflashCt_choice = zeros(ns,nt_choice);
   rflashCt_end = zeros(ns,nt_end);
   rflashCt_pokeend = zeros(ns,nt_pokeend);
   
   lbupsCt_start = zeros(ns,nt_start);
   lbupsCt_lastflash = zeros(ns,nt_lastflash);
   lbupsCt_choice = zeros(ns,nt_choice);
   lbupsCt_end = zeros(ns,nt_end);
   lbupsCt_pokeend = zeros(ns,nt_pokeend);
   
   rbupsCt_start = zeros(ns,nt_start);
   rbupsCt_lastflash = zeros(ns,nt_lastflash);
   rbupsCt_choice = zeros(ns,nt_choice);
   rbupsCt_end = zeros(ns,nt_end);
   rbupsCt_pokeend = zeros(ns,nt_pokeend);
   
   for k = 1:ns
       
      %bin lflash relative to event
      if numel(handles.lflashes{k}) > 0
          lflashCt_start(k,:) = sum(discretize(handles.lflashes{k},xvec_start)== (1:nt_start)',2);
          lflashCt_lastflash(k,:) = sum(discretize(handles.lflashes{k},xvec_lastflash)== (1:nt_lastflash)',2);
          lflashCt_choice(k,:) = sum(discretize(handles.lflashes{k},xvec_choice)== (1:nt_choice)',2);
          lflashCt_end(k,:) = sum(discretize(handles.lflashes{k},xvec_end)== (1:nt_end)',2);
          lflashCt_pokeend(k,:) = sum(discretize(handles.lflashes{k},xvec_pokeend)== (1:nt_pokeend)',2);
      end
      
      %bin rflash relative to event
      if numel(handles.rflashes{k}) > 0
          rflashCt_start(k,:) = sum(discretize(handles.rflashes{k},xvec_start)== (1:nt_start)',2);
          rflashCt_lastflash(k,:) = sum(discretize(handles.rflashes{k},xvec_lastflash)== (1:nt_lastflash)',2);
          rflashCt_choice(k,:) = sum(discretize(handles.rflashes{k},xvec_choice)== (1:nt_choice)',2);
          rflashCt_end(k,:) = sum(discretize(handles.rflashes{k},xvec_end)== (1:nt_end)',2);
          rflashCt_pokeend(k,:) = sum(discretize(handles.rflashes{k},xvec_pokeend)== (1:nt_pokeend)',2);
      end

      %bin lbups relative to event
      if numel(handles.lbups{k}) > 0
          lbupsCt_start(k,:) = sum(discretize(handles.lbups{k},xvec_start)== (1:nt_start)',2);
          lbupsCt_lastflash(k,:) = sum(discretize(handles.lbups{k},xvec_lastflash)== (1:nt_lastflash)',2);
          lbupsCt_choice(k,:) = sum(discretize(handles.lbups{k},xvec_choice)== (1:nt_choice)',2);
          lbupsCt_end(k,:) = sum(discretize(handles.lbups{k},xvec_end)== (1:nt_end)',2);
          lbupsCt_pokeend(k,:) = sum(discretize(handles.lbups{k},xvec_pokeend)== (1:nt_pokeend)',2);
      end
      
      %bin rbups relative to event
      if numel(handles.rbups{k}) > 0
          rbupsCt_start(k,:) = sum(discretize(handles.rbups{k},xvec_start)== (1:nt_start)',2);
          rbupsCt_lastflash(k,:) = sum(discretize(handles.rbups{k},xvec_lastflash)== (1:nt_lastflash)',2);
          rbupsCt_choice(k,:) = sum(discretize(handles.rbups{k},xvec_choice)== (1:nt_choice)',2);
          rbupsCt_end(k,:) = sum(discretize(handles.rbups{k},xvec_end)== (1:nt_end)',2);
          rbupsCt_pokeend(k,:) = sum(discretize(handles.rbups{k},xvec_pokeend)== (1:nt_pokeend)',2);
      end      
      
   end
   
   stim = struct();
   stim.lflashCt_start = lflashCt_start;
   stim.lflashCt_lastflash = lflashCt_lastflash;
   stim.lflashCt_choice = lflashCt_choice;
   stim.lflashCt_end = lflashCt_end;
   stim.lflashCt_pokeend = lflashCt_pokeend;
   
   stim.rflashCt_start = rflashCt_start;
   stim.rflashCt_lastflash = rflashCt_lastflash;
   stim.rflashCt_choice = rflashCt_choice;
   stim.rflashCt_end = rflashCt_end;
   stim.rflashCt_pokeend = rflashCt_pokeend;
   
   stim.lbupsCt_start = lbupsCt_start;
   stim.lbupsCt_lastflash = lbupsCt_lastflash;
   stim.lbupsCt_choice = lbupsCt_choice;
   stim.lbupsCt_end = lbupsCt_end;
   stim.lbupsCt_pokeend = lbupsCt_pokeend;
   
   stim.rbupsCt_start = rbupsCt_start;
   stim.rbupsCt_lastflash = rbupsCt_lastflash;
   stim.rbupsCt_choice = rbupsCt_choice;
   stim.rbupsCt_end = rbupsCt_end;
   stim.rbupsCt_pokeend = rbupsCt_pokeend;
   
   A{j}.stim = stim;
   %}



