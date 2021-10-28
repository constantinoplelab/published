%script to scrub and simplify the ofc data

%will remove the following neurons via concatdata_ofc_all:
%   multiunit
%   supremely low firing rate (failed t-test for usability)

%this script will do the following
% remov low firing rate < 1Hz
% will remove the following opt-out trials
% will omit the spike times

% will include the following additional details:
%   a trial index to show which trials were kept in the scrubbed data

%load concatenated form
savedir = '~/projects/ofc/data/published/'; 
%savedir = '~/projects/ofc/data/processed/'; %for my data

%savedir = '~/projects/ofc/data/PhysiologyData/Christine_OFC_2019';
fname_load = strcat(savedir,'concatdata_ofc.mat'); %saved file
fname_save = strcat(savedir,'concatdata_ofc_noOptOut.mat');
B = load(fname_load);
A = B.A; %easier to rename structure here

nf = numel(A);

keepvec = zeros(nf,1); %for removing low firing rate neurons

for j = 1:nf
    
    if mod(j,100) == 0
        disp(j)
    end
    
    if mean(mean(A{j}.hmat_start)) < 1 %average firing rate less than one Hz, then remove
        continue
    else
       keepvec(j) = 1;

       hits = A{j}.hits;
       mask = ~isnan(hits); %which files to keep

       trial_idx = find(mask);

       A{j}.hmat_start = A{j}.hmat_start(mask,:);
       A{j}.hmat_pokeend = A{j}.hmat_pokeend(mask,:);
       A{j}.hmat_choice = A{j}.hmat_choice(mask,:);
       A{j}.hmat_lastflash = A{j}.hmat_lastflash(mask,:);
       A{j}.hmat_end = A{j}.hmat_end(mask,:);
       A{j}.hmat_leavecpoke = A{j}.hmat_leavecpoke(mask,:);

       A{j}.handles.start = A{j}.handles.start(mask);
       A{j}.handles.end = A{j}.handles.end(mask);
       A{j}.handles.choice = A{j}.handles.choice(mask);
       A{j}.handles.leavecpoke = A{j}.handles.leavecpoke(mask);  

       A{j}.hits = hits(mask);
       A{j}.prev_hits = A{j}.prev_hits(mask);
       A{j}.went_right = A{j}.went_right(mask); %different than parse_choice
       A{j}.this_left_volume = A{j}.this_left_volume(mask);
       A{j}.this_right_volume = A{j}.this_right_volume(mask);
       A{j}.left_prob = A{j}.left_prob(mask);
       A{j}.right_prob = A{j}.right_prob(mask);
       A{j}.chosenprob = A{j}.chosenprob(mask);
       A{j}.chosenval = A{j}.chosenval(mask);

       A{j}.nspikes = A{j}.nspikes(mask);

       A{j}.trial_idx = trial_idx;

    end
    
end

xvec_start = B.xvec_start;
xvec_pokeend = B.xvec_pokeend;
xvec_choice = B.xvec_choice;
xvec_lastflash = B.xvec_lastflash;
xvec_end = B.xvec_end;
xvec_leavecpoke = B.xvec_leavecpoke;

A = A(find(keepvec));


save(fname_save,'A','xvec*','-v7.3');

%% if flash and beep timing included

%if flashes and beeps timing added in
%A{j}.handles.lflashes =  A{j}.handles.lflashes(mask);
%A{j}.handles.rflashes =  A{j}.handles.rflashes(mask);
%A{j}.handles.lbups =  A{j}.handles.lbups(mask);
%A{j}.handles.rbups =  A{j}.handles.rbups(mask);
   
%handle the stim
   %{
   A{j}.stim.lflashCt_start = A{j}.stim.lflashCt_start(mask);
   A{j}.stim.lflashCt_lastflash = A{j}.stim.lflashCt_lastflash(mask);
   A{j}.stim.lflashCt_choice = A{j}.stim.lflashCt_choice(mask);
   A{j}.stim.lflashCt_end = A{j}.stim.lflashCt_end(mask);
   A{j}.stim.lflashCt_pokeend = A{j}.stim.lflashCt_pokeend(mask);
   
   A{j}.stim.rflashCt_start = A{j}.stim.rflashCt_start(mask);
   A{j}.stim.rflashCt_lastflash = A{j}.stim.rflashCt_lastflash(mask);
   A{j}.stim.rflashCt_choice = A{j}.stim.rflashCt_choice(mask);
   A{j}.stim.rflashCt_end = A{j}.stim.rflashCt_end(mask);
   A{j}.stim.rflashCt_pokeend = A{j}.stim.rflashCt_pokeend(mask);
   
   A{j}.stim.lbupsCt_start = A{j}.stim.lbupsCt_start(mask);
   A{j}.stim.lbupsCt_lastflash = A{j}.stim.lbupsCt_lastflash(mask);
   A{j}.stim.lbupsCt_choice = A{j}.stim.lbupsCt_choice(mask);
   A{j}.stim.lbupsCt_end = lA{j}.stim.lbupsCt_end(mask);
   A{j}.stim.lbupsCt_pokeend = A{j}.stim.lbupsCt_pokeend(mask);
   
   A{j}.stim.rbupsCt_start = A{j}.stim.rbupsCt_start(mask);
   A{j}.stim.rbupsCt_lastflash = A{j}.stim.rbupsCt_lastflash(mask);
   A{j}.stim.rbupsCt_choice = A{j}.stim.rbupsCt_choice(mask);
   A{j}.stim.rbupsCt_end = A{j}.stim.rbupsCt_end(mask);
   A{j}.stim.rbupsCt_pokeend = A{j}.stim.rbupsCt_pokeend(mask);
   %}
d