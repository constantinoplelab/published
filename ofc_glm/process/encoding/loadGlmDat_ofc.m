function [Dat,stim] = loadGlmDat_ofc(ops,datname)
%parses ofc data into Dat structure for glm optimization
%
%for ns trials and nt data points in trial-based, or ns*nt data samples
%with ncov covariates, Dat has required fields:
%
%Dat.y: (1,ns*nt) or (nt,ns) matrix of spike counts
%Dat.trial_idx: (1,ns) set of trial numbers
%Dat.start_idx = (1,ns) vector of indices for start of trials
%Dat.end_idx = (1,ns) vector of indices for end of trials
%Dat.covariates = (1,ncov-1) cell array of covariates
%Dat.covnames = (1,ncov-1) cell array of covariate names
%
%Stim struct will have the following fields:
%
%stim.TDstim. (1,nc) cell of stimuli that evoke a kernel response
%stim.TDcausalvec. (2,nc) matrix of booleans for causal or acausal kernel
%response
%stim.xmat. (1,ncov-nc) cell of time-independent covariates. bias terms
%stim.Xorder. labels for stim.xmat
%stim.stimlegend. labels for TDstim

%% load glm data

dispL('loading data',5,ops.dispL);
ofcdat = load(datname);

%trial hits, chosen prob and chosen val
hits = ofcdat.hits;
prev_hits = [nan;hits(1:end-1)]; 
ofcdat.prev_hits = prev_hits;
ofcdat.prev_optouts = isnan(prev_hits); %previous trial was opt-out

%chosen volume

chosenval = (ofcdat.went_right==1).*ofcdat.this_right_volume + ...
    (ofcdat.went_right==0).*ofcdat.this_right_volume;
chosenval(isnan(ofcdat.went_right)) = nan;


%scrub data for duplicate values and nans. remove those trials

%trial start, leave cpoke, choice. remove from trial_idx
%check rogue nan values.
kmask1 =~isnan(ofcdat.handles.start);
kmask2 = ~isnan(ofcdat.handles.leavecpoke);
kmask3 = ~isnan(ofcdat.handles.choice);

%check for duplicate values. add to mask to filter out duplicates
[~,u1] = unique(ofcdat.handles.start);
[~,u2] = unique(ofcdat.handles.leavecpoke);
[~,u3] = unique(ofcdat.handles.choice);
kmask1a = false(size(ofcdat.handles.start));
kmask1a(u1) = true;
kmask2a = false(size(ofcdat.handles.leavecpoke));
kmask2a(u2) = true;
kmask3a = false(size(ofcdat.handles.choice));
kmask3a(u3) = true;

kmask1 = kmask1 & kmask1a;
kmask2 = kmask2 & kmask2a;
kmask3 = kmask3 & kmask3a;

%scrub any bad data (duplicate value, nans) from all of the raw data, but
%leave opt-out trials. might be needed later, and it is not technically
%"bad" data, just unused for this task
kmask = kmask1 & kmask2 & kmask3;

ofcdat.handles.start = ofcdat.handles.start(kmask);
ofcdat.handles.leavecpoke = ofcdat.handles.leavecpoke(kmask);
ofcdat.handles.choice = ofcdat.handles.choice(kmask);
ofcdat.handles.lflashes = ofcdat.handles.lflashes(kmask);
ofcdat.handles.rflashes = ofcdat.handles.rflashes(kmask);
ofcdat.handles.lbups = ofcdat.handles.lbups(kmask);
ofcdat.handles.rbups = ofcdat.handles.rbups(kmask);
ofcdat.hits = ofcdat.hits(kmask);
ofcdat.prev_hits = ofcdat.prev_hits(kmask);
ofcdat.prev_optouts = ofcdat.prev_optouts(kmask);
ofcdat.right_prob = ofcdat.right_prob(kmask);
ofcdat.left_prob = ofcdat.left_prob(kmask);
ofcdat.went_right = ofcdat.went_right(kmask);

chosenval = chosenval(kmask);


%removing the opt-out trials from list of usable trials 
if isfield(ops,'removeData')
    if ops.removeData
        dispL('removing opt out trials',5,ops.dispL)
        mask = ~isnan(ofcdat.hits);       
    else
        mask = true(size(ofcdat.hits));       
    end  
end

%make trial vector to save for later. omits opt-out trials
trial_idx = 1:numel(ofcdat.handles.start);
trial_idx = trial_idx(mask);

%% parse spikes and stimuli into bins

%find smallest and largest event time
tmin = min([ofcdat.spiketimes',[ofcdat.handles.lflashes{:}],[ofcdat.handles.rflashes{:}],...
    [ofcdat.handles.lbups{:}],[ofcdat.handles.rbups{:}],ofcdat.handles.start',...
    ofcdat.handles.leavecpoke',ofcdat.handles.choice']);
tmax = max([ofcdat.spiketimes',[ofcdat.handles.lflashes{:}],[ofcdat.handles.rflashes{:}],...
    [ofcdat.handles.lbups{:}],[ofcdat.handles.rbups{:}],ofcdat.handles.start',...
    ofcdat.handles.leavecpoke',ofcdat.handles.choice']);

tmesh = tmin:ops.dt:tmax+ops.dt;
nt = numel(tmesh)-1;

%spike times
y = histcounts(ofcdat.spiketimes,tmesh);

%flashes and beeps
lflash_vec = [ofcdat.handles.lflashes{:}]; 
rflash_vec = [ofcdat.handles.rflashes{:}]; 
lbups_vec = [ofcdat.handles.lbups{:}]; 
rbups_vec = [ofcdat.handles.rbups{:}]; 

stim_lflash = histcounts(lflash_vec,tmesh)'; 
stim_rflash = histcounts(rflash_vec,tmesh)'; 
stim_lbups = histcounts(lbups_vec,tmesh)'; 
stim_rbups = histcounts(rbups_vec,tmesh)';

stim_start = histcounts(ofcdat.handles.start,tmesh);
stim_leavecpoke = histcounts(ofcdat.handles.leavecpoke,tmesh);
stim_choice = histcounts(ofcdat.handles.choice,tmesh);

cpoke_idx = find(stim_leavecpoke);
choice_idx = find(stim_choice);


%decide on trial start and end. currently 2s before enter center poke until
%2s before next trial start
start_idx = max([1,choice_idx(1:end-1)+1],find(stim_start)-ceil(2.0/ops.dt)); %try allowing a 2s lead
end_idx = min(ceil(tmax/ops.dt),[start_idx(2:end)-1,nt]); %trial should endure right up until next trial starts

%start of the trial w.r.t. animal and covariates, when animal enters cpoke
entercpoke_idx = find(stim_start);


%% parameters of interest

%trial-level parameters
hits = ofcdat.hits;
prev_hits = ofcdat.prev_hits;

%binary behaviors 
isRightSafe=  ofcdat.right_prob==1;
isLeftSafe=   ofcdat.left_prob==1;   

choseSafe = (isRightSafe & ofcdat.went_right==1) | ...
        (isLeftSafe & ofcdat.went_right==0);
XprevReward = prev_hits ==1  ;
XprevLoss = prev_hits ==0 ;   
Xreward = hits == 1;
Xloss =   hits == 0;
Xsafe = choseSafe & ~isnan(hits); %consider change to hits==1.
Xrisky = ~choseSafe & ~isnan(hits);
Xleft = ofcdat.went_right==0;
Xright = ofcdat.went_right==1;

covbias = ones(nt,1); %background firing rate

%rewarded water and expected reward
water = hits.*chosenval(1:length(hits)); %water received on each trial
water(isnan(water))=0; %omit nans from opt outs. change to no water
water_prev = [0;water(1:end-1)]; %volume rewarded on previous trial

%average expected reward as running average, the rewardRate  
rewRate = cumsum(water,'omitnan')./(1:numel(water))'; 
rewRate(isnan(rewRate))=0; %some nans at the beginning
prevRewRate = [0;rewRate(1:end-1)]; %average reward up until previous trial

%session progress
session_progress = (1:numel(ofcdat.handles.start))/(numel(ofcdat.handles.start));

%% parameters converted to covariates by being time-locked to event

%combine trial-level info with event stimuli to make covariates
covPrevWin = zeros(nt,1); 
covPrevWin(entercpoke_idx(XprevReward)) = 1;

covPrevLoss = zeros(nt,1);
covPrevLoss(entercpoke_idx(XprevLoss)) = 1;

covLeft = zeros(nt,1);
covLeft(cpoke_idx(Xleft)) = 1;

covRight = zeros(nt,1);
covRight(cpoke_idx(Xright)) = 1;

covWin = zeros(nt,1);
covWin(choice_idx(Xreward)) = 1;

covLoss = zeros(nt,1);
covLoss(choice_idx(Xloss)) = 1;

covRisky = zeros(nt,1);
covRisky(choice_idx(Xrisky)) = 1;

covSafe = zeros(nt,1);
covSafe(choice_idx(Xsafe)) = 1;

%new ones. ---------------
%value of previous water instead of previous win/loss
covPrevVol = zeros(nt,1);
covPrevVol(entercpoke_idx(XprevReward)) = water_prev(XprevReward)/48; %scale by max value

%water on current trial, excluding 0 ul
covVol_no0 = zeros(nt,1);
covVol_no0(choice_idx(Xreward)) = water(Xreward)/48;

%prevReward rate
covPrevRewardRate = zeros(nt,1);
covPrevRewardRate(entercpoke_idx) = prevRewRate/48;

%a previous opt-out kernel
covPrevOptOut= zeros(nt,1);
covPrevOptOut(entercpoke_idx(ofcdat.prev_optouts)) = 1; %scale by max value

%session progress as proxy for motivation
covSessProg = zeros(nt,1);
covSessProg(entercpoke_idx) = session_progress;


%% make covariates based on stimtype

switch ops.stimtype
    
    
    case '20200710'   
        %all stimuli are fully time-dependent,to allow for complete overlap
        %from events.
        %
        %start prev. win/loss
        %(acausal) start prev. win/loss
        %leave cpoke left/right
        %choice (anticausal) safe/risky
        %choice (causal) win lose
        %time-dep. stimuli parametrized as left/right clicks and flashes
        
        covariates = {stim_lflash,stim_rflash,stim_lbups, stim_rbups,...
            covPrevWin, covPrevLoss,covLeft, covRight,...
            covWin, covLoss,covRisky, covSafe,...
            covbias};

        
        TDcausalvec = [1,1,1,1,1,1,1,1,1,1,0,0; %row1: causal kernels
                       0,0,0,0,1,1,0,0,0,0,1,1]; %row2: acausal kernel
                   
        
        covnames = {'LFlash','RFlash','LClick','RClick',...
        'prevWin','prevLoss','Left','Right',...
        'win','loss','prevWin (ac)','prevLoss (ac)',...
        'risky (ac)','safe (ac)'};

    
    case '20200903'   
        % "best" model used in paper.

        covariates = {stim_lflash,stim_rflash,stim_lbups, stim_rbups,...
            covPrevWin, covPrevLoss,covLeft, covRight,...
            covWin, covLoss,covRisky, covSafe,...
            covPrevRewardRate,covPrevOptOut,covSessProg,...
            covbias};

        
        TDcausalvec = [1,1,1,1,1,1,1,1,1,1,0,0,1,1,1; %row1: causal filters
                       0,0,0,0,1,1,0,0,0,0,1,1,1,1,1]; %row2: anticausal choice as well 
                   
        
        covnames = {'LFlash','RFlash','LClick','RClick',...
        'prevWin','prevLoss','Left','Right',...
        'win','loss',...
        'prevRewardRate',...
        'prevOptOut','sessProg'... %end of causal
        'prevWin (ac)','prevLoss (ac)',...     
        'risky (ac)','safe (ac)',...
        'prevRewardRate (ac)','prevOptOut (ac)','sessProg (ac)'};  
    
    case '20200930-1'   
        % "+prev_vol" in paper
        % 2020903 + prev_vol - prevWin

        covariates = {stim_lflash,stim_rflash,stim_lbups, stim_rbups,...
            covPrevLoss,covLeft, covRight,...
            covWin, covLoss,covRisky, covSafe,...
            covPrevRewardRate,covPrevOptOut,covSessProg,covPrevVol,...
            covbias};

        
        TDcausalvec = [1,1,1,1,1,1,1,1,1,0,0,1,1,1,1; %row1: causal filters
                       0,0,0,0,1,0,0,0,0,1,1,1,1,1,1]; %row2: anticausal choice as well 
                   
        
        covnames = {'LFlash','RFlash','LClick','RClick',...
        'prevLoss','Left','Right',...
        'win','loss',...
        'prevRewardRate',...
        'prevOptOut','sessProg',...
        'prevVol',...%end of causal
        'prevLoss (ac)',...     
        'risky (ac)','safe (ac)',...
        'prevRewardRate (ac)','prevOptOut (ac)','sessProg (ac)','prevVol (ac)'};
    
    case '20200930-2'   
        %"+vol" in paper
        % 2020903 + prev_vol - prevWin + vol_no0 - win

        covariates = {stim_lflash,stim_rflash,stim_lbups, stim_rbups,...
            covPrevLoss,covLeft, covRight,...
            covLoss,covRisky, covSafe,...
            covPrevRewardRate,covPrevOptOut,covSessProg,covPrevVol,covVol_no0...
            covbias};

        
        TDcausalvec = [1,1,1,1,1,1,1,1,0,0,1,1,1,1,1; %row1: causal filters
                       0,0,0,0,1,0,0,0,1,1,1,1,1,1,0]; %row2: anticausal choice as well 
                   
        
        covnames = {'LFlash','RFlash','LClick','RClick',...
        'prevLoss','Left','Right',...
        'loss',...
        'prevRewardRate',...
        'prevOptOut','sessProg',...
        'prevVol','vol_no0'...%end of causal
        'prevLoss (ac)',...     
        'risky (ac)','safe (ac)',...
        'prevRewardRate (ac)','prevOptOut (ac)','sessProg (ac)','prevVol (ac)'};

    
    case '20210302'   
        %"-rewHist" in paper
        %20200710 - prevWin - prevLoss
        
        covariates = {stim_lflash,stim_rflash,stim_lbups, stim_rbups,...
            covLeft, covRight,...
            covWin, covLoss,covRisky, covSafe,...
            covbias};

        
        TDcausalvec = [1,1,1,1,1,1,1,1,0,0; %row1: causal kernels
                       0,0,0,0,0,0,0,0,1,1]; %row2: acausal kernel
                   
        
        covnames = {'LFlash','RFlash','LClick','RClick',...
        'Left','Right',...
        'win','loss',...
        'risky (ac)','safe (ac)'};
    
    otherwise
        disp('you have not entered a known stimulus configuration')
end

nkern = sum(sum(TDcausalvec));

Dat = struct();
Dat.y = y;
Dat.start_idx = start_idx;
Dat.end_idx = end_idx;
Dat.covariates = covariates;
Dat.covnames = covnames;
Dat.ofcdat = ofcdat;
%specific alignment stuff for ofc
Dat.entercpoke_idx = entercpoke_idx;
Dat.choice_idx = choice_idx;
Dat.leavecpoke_idx = cpoke_idx;
Dat.nkern = nkern;

if strcmp(ops.stimtype,'20200819-3')
    Dat.V = V;
    Dat.covnames_orig = covnames_orig;
end

%check if ordering of events is correct in each trial. if not, omit trial
omitmask = false(size(trial_idx));
for j = 1:numel(trial_idx)
    tj = trial_idx(j);
    check = diff([start_idx(tj), entercpoke_idx(tj), cpoke_idx(tj), choice_idx(tj), end_idx(tj)]);
    if sum(check < 0) > 0
        omitmask(j) = true;
    end
        
end
Dat.trial_idx = trial_idx(~omitmask);

%stimulus data structure. Not sure if this is redundant
stim = struct();      
stim.TDstim = Dat.covariates(1:end-1); %time-dependent covariates       
stim.TDcausalvec = TDcausalvec;  
stim.xmat = Dat.covariates(end); %time-independent covariates. but still stored as time-dep. for ease of LL cal.
%label information
stim.Xorder = {'bias'};
stim.stimlegend = Dat.covnames;
end

