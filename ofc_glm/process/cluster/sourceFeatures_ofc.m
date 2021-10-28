function [fname] = sourceFeatures_ofc(ops)
    %generates building blocks of features from OFC data.
    
    
    %options for how the feature coding will be done
    nsamp = ops.nsamps; %how many bootstrap for ARI samples will be made?
    useZscore = ops.usezscore;
    bootstrap_proportion = ops.bootstrap_proportion; %proportion of samples in bootstrap
    bootstrap_replacemet = ops.bootstrap_replacement; % with replacement= = true
    normfet = ops.normfet; %is z-score is false, should features be normalzied by max response?
    datadir = ops.datadir; %where will data be stored
    
    %how will this data combo be notated    
    if bootstrap_replacemet
        bsrep = 'Rep_';
    else
        bsrep = 'noRep_';
    end
    
    if useZscore
        ztext = 'Z_';
    else
        if normfet
            ztext = 'Norm_';
        else
            ztext = 'noNorm_';
        end
    end
    
    nstext = strcat('ns_',num2str(nsamp));
    
    bstext = strcat('bs',num2str(round(bootstrap_proportion*100)));
    
    fname = strcat('ofcFeaturebase_',ztext,bstext,bsrep,nstext,'.mat');
    
    
    
%% load cleaned data

disp('loading concatenated ephys data')

D = load(fullfile(ops.datadir,'/concatdata_ofc_noOptOut.mat'));
A = D.A; %the main data structure for saved data
nn = numel(A);

Aj = 1:nn; %relic indexer of largeer dataset
        
%% initialize features

%forms of presented value offers and rewarded values
volvec = [6,12,24,48]; %offered volume, via clicks
probvec = 0:0.1:1; %offered probability, via flashes
EV_bins = [0,6,12,24,48]; %prob * vol
EV_bins_coarse = [0,0.01,12,48]; %0,low,high
RPE_bins = [-fliplr(EV_bins(2:end)),0,EV_bins(2:end)]; %RPE of reward
DEV_bins = RPE_bins; %EV_left - EV_right
%RPE_bins = [-fliplr(EV_bins(3:end)),0,EV_bins(2:end)];
%EV_all = unique(volvec'*probvec);
%RPE_all = unique(volvec-EV_all);
%EV_bins = 2.^(linspace(log2(0.6),log2(48),8));

%size info about base value vectors
nCO = 8; %number of features in marginal win/loss, prevWin/Loss, L/R, riskySafe
nbase_vol = numel(volvec);
nbase_prob = numel(probvec); %0:0.1:1
nbase_EV = numel(EV_bins)-1;
nbase_EV_coarse = numel(EV_bins_coarse)-1;
nbase_RPE = numel(RPE_bins)-1;
nbase_DEV = numel(RPE_bins)-1;
nbase_rew = numel(EV_bins);
nbase_rew_coarse = numel(EV_bins_coarse);

%choice and outcome features
%TODO: remove
COfet = cell(nCO,1);
COfet_sem = cell(nCO,1);
COfet_legend = {'win','loss','prevWin','prevLoss','left','right','risky','safe'};
for j = 1:nCO
    COfet{j} = zeros(nn,1,nsamp+1);
    COfet_sem{j} = zeros(nn,1);
end

%reward offer features
volvecChoose = zeros(nn,nbase_vol,nsamp+1);
probvecChoose = zeros(nn,nbase_prob,nsamp+1);
EVvecChoose = zeros(nn,nbase_EV,nsamp+1);
RPEvec = zeros(nn,nbase_RPE,nsamp+1);
rewVolvec = zeros(nn,nbase_rew,nsamp+1);

volvecChoose_sem = zeros(nn,nbase_vol); %3. all data and  two test sets
probvecChoose_sem = zeros(nn,nbase_prob);
EVvecChoose_sem = zeros(nn,nbase_EV);
RPEvec_sem = zeros(nn,nbase_RPE);
rewVolvec_sem = zeros(nn,nbase_rew);


%presented offer features
volvecL = zeros(nn,nbase_vol,nsamp+1);
probvecL = zeros(nn,nbase_prob,nsamp+1);
EVvecL = zeros(nn,nbase_EV,nsamp+1);
volvecR = zeros(nn,nbase_vol,nsamp+1); 
probvecR = zeros(nn,nbase_prob,nsamp+1);
EVvecR = zeros(nn,nbase_EV,nsamp+1);
DEVvec = zeros(nn,nbase_DEV,nsamp+1);

volvecL_sem = zeros(nn,nbase_vol);
probvecL_sem = zeros(nn,nbase_prob);
EVvecL_sem = zeros(nn,nbase_EV);
volvecR_sem = zeros(nn,nbase_vol); 
probvecR_sem = zeros(nn,nbase_prob);
EVvecR_sem = zeros(nn,nbase_EV);
DEVvec_sem = zeros(nn,nbase_DEV);

%safe/risky parametriztion of presented value offers
volvecSafe = zeros(nn,nbase_vol,nsamp+1); 
probvecSafe = zeros(nn,nbase_prob,nsamp+1);
EVvecSafe = zeros(nn,nbase_EV,nsamp+1);
volvecRisky = zeros(nn,nbase_vol,nsamp+1);
probvecRisky = zeros(nn,nbase_prob,nsamp+1);
EVvecRisky = zeros(nn,nbase_EV,nsamp+1);

volvecSafe_sem = zeros(nn,nbase_vol); 
probvecSafe_sem = zeros(nn,nbase_prob);
EVvecSafe_sem = zeros(nn,nbase_EV);
volvecRisky_sem = zeros(nn,nbase_vol); 
probvecRisky_sem = zeros(nn,nbase_prob);
EVvecRisky_sem = zeros(nn,nbase_EV);

%% set trials mask and time window masks for averaging each feature


%create function handles for the sub masks
mask_cell = cell(nCO);
mask_cell{1} = @(x) A{x}.hits == 1; %win
mask_cell{2} = @(x) A{x}.hits == 0; %loss
mask_cell{3} = @(x) A{x}.prev_hits == 1; %previous win
mask_cell{4} = @(x) A{x}.prev_hits == 0; %previous loss
mask_cell{5} = @(x) A{x}.went_right==0; %left choice
mask_cell{6} = @(x) A{x}.went_right==1; %right choice
mask_cell{7} = @(x) A{x}.chosenprob~=1; %risky choice
mask_cell{8} = @(x) A{x}.chosenprob==1; %safe choice

rng(101)

%set the time window and alignment to average firing rate for each type of feature
align_reward = 'choice';
align_PO = 'choice';
align_CO = {'choice','choice','start','start',...
    'leavecpoke','leavecpoke','leavecpoke','leavecpoke'};
wind_reward = [0,3]; %rewarded offer features
%wind_PO = [-4,0];
wind_PO = [-1,0];
wind_CO = {[0,3],[0,3],[-1,2],[-1,2],[0,1.5],[0,1.5],[0,1.5],[0,1.5]};

%create and save a source options to recall in plotting
sourceops = struct();
sourceops.align_reward = align_reward;
sourceops.align_PO = align_PO;
sourceops.align_CO = align_CO;
sourceops.wind_reward = wind_reward;
sourceops.wind_PO = wind_PO;
sourceops.wind_CO = wind_CO;
sourceops.volvec = volvec;
sourceops.probvec = probvec;
sourceops.EV_bins = EV_bins;
sourceops.EV_bins_coarse = EV_bins_coarse;
sourceops.RPE_bins = RPE_bins;
sourceops.DEV_bins = DEV_bins;




%% build features

disp('building features')
for jj = 1:nn
    
    j = Aj(jj); %the actual index
    
    %disp(j)

    %get max firing rate of cell to normalize
    if normfet
        nf = nanmax(A{j}.hmat_start,[],'all');
    else
        nf = 1;
    end
    meandat = nanmean(nanmean(A{j}.hmat_start)); %unsure if I should use mean of mean, or mean of all?
    stddat = nanstd(nanmean(A{j}.hmat_start));

    %bootstrap masks
    ns = numel(A{j}.nspikes); %numbers total samples
    tmat = false(ns,nsamp);
    for m = 1:nsamp
        randoff = randsample(ns,floor(bootstrap_proportion*ns),bootstrap_replacemet);
        tmat(randoff,m) = true;
    end

    %choose heatmats for reward and presented offer features
    switch align_reward
        case 'start'
            hmat_reward = A{j}.hmat_start;
            tmesh_reward = -2:0.05:4;
        case 'leavecpoke'
            hmat_reward = A{j}.hmat_leavecpoke;
            tmesh_reward = -4:0.05:4;
        case 'choice'
            hmat_reward = A{j}.hmat_choice;     
            tmesh_reward = -4:0.05:4;
    end       
    %grab time windows of interest and average
    wind_idx = discretize(wind_reward,tmesh_reward);
    hmat_reward = nanmean(hmat_reward(:,wind_idx(1):wind_idx(2)),2);

    %choose heatmats for reward and presented offer features
    switch align_PO
        case 'start'
            hmat_PO = A{j}.hmat_start;
            tmesh_PO = -2:0.05:4;
        case 'leavecpoke'
            hmat_PO = A{j}.hmat_leavecpoke;
            tmesh_PO = -4:0.05:4;
        case 'choice'
            hmat_PO = A{j}.hmat_choice;     
            tmesh_PO = -4:0.05:4;
    end       
    %grab time windows of interest and average
    wind_idx = discretize(wind_PO,tmesh_PO);
    hmat_PO = nanmean(hmat_PO(:,wind_idx(1):wind_idx(2)),2);


    %build the choice outcome features
    for m = 1:nCO
        mask =  mask_cell{m}(j);

        %decide which heatmat to use
        switch align_CO{m}
            case 'start'
                hmat = A{j}.hmat_start;
                tmesh = -2:0.05:4;
            case 'leavecpoke'
                hmat = A{j}.hmat_leavecpoke;
                tmesh = -4:0.05:4;
            case 'choice'
                hmat = A{j}.hmat_choice;     
                tmesh = -4:0.05:4;
        end

        %average firing rate in that time window, all trials
        wind_idx = discretize(wind_CO{m},tmesh);
        hmat = nanmean(hmat(:,wind_idx(1):wind_idx(2)),2);

        if useZscore
            COfet{m}(jj,1) = nanmean(hmat(mask)-meandat)/stddat;
            COfet_sem{m}(jj) = nanstd((hmat(mask)-meandat)/stddat)/sqrt(numel(hmat(mask)));

            for l = 1:nsamp
                COfet{m}(jj,l+1) = nanmean(hmat(mask & tmat(:,l))-meandat)/stddat;
            end
        else           
            COfet{m}(jj,1) = nanmean(hmat(mask))/nf;
            COfet_sem{m}(jj) = nanstd(hmat(mask)/nf)/sqrt(numel(hmat(mask)));

            for l = 1:nsamp
                 COfet{m}(jj,l+1) = nanmean(hmat(mask & tmat(:,l)))/nf;
            end

        end
    end


    %build the volume features

    %reparametrize left right into into safe and risky volumes
    safeOptL = A{j}.left_prob ==1 ;
    safeOptR = A{j}.right_prob==1;
    riskyOptL = A{j}.left_prob ~=1 ;
    riskyOptR = A{j}.right_prob ~=1;

    for k = 1:nbase_vol      
        mask = A{j}.chosenval==volvec(k); %chosen value
        maskL = A{j}.this_left_volume==volvec(k);
        maskR = A{j}.this_right_volume==volvec(k);
        maskSafe =  (A{j}.this_left_volume==volvec(k) & safeOptL) | ...
                    (A{j}.this_right_volume==volvec(k) & safeOptR);
        maskRisky = (A{j}.this_left_volume==volvec(k) & riskyOptL) | ...
                    (A{j}.this_right_volume==volvec(k) & riskyOptR);    

        if useZscore
            volvecChoose(jj,k,1) = nanmean(hmat_reward(mask)-meandat)/stddat;
            volvecL(jj,k,1) = nanmean(hmat_PO(maskL)-meandat)/stddat;
            volvecR(jj,k,1) = nanmean(hmat_PO(maskR)-meandat)/stddat;
            %use left-right hmats. they are what we want
            volvecSafe(jj,k,1) = nanmean(hmat_PO(maskSafe)-meandat)/stddat;
            volvecRisky(jj,k,1) = nanmean(hmat_PO(maskRisky)-meandat)/stddat;
            %sem
            volvecChoose_sem(jj,k) = nanstd((hmat_reward(mask)-meandat)/stddat)/sqrt(sum(mask));
            volvecL_sem(jj,k) = nanstd((hmat_PO(maskL)-meandat)/stddat)/sqrt(sum(maskL));
            volvecR_sem(jj,k) = nanstd((hmat_PO(maskR)-meandat)/stddat)/sqrt(sum(maskR));
            volvecSafe_sem(jj,k) = nanstd((hmat_PO(maskSafe)-meandat)/stddat)/sqrt(sum(maskSafe));
            volvecRisky_sem(jj,k) = nanstd((hmat_PO(maskRisky)-meandat)/stddat)/sqrt(sum(maskRisky));

            for l = 1:nsamp
                volvecChoose(jj,k,l+1) = nanmean(hmat_reward(mask & tmat(:,l))-meandat)/stddat;
                volvecL(jj,k,l+1) = nanmean(hmat_PO(maskL & tmat(:,l))-meandat)/stddat;
                volvecR(jj,k,l+1) = nanmean(hmat_PO(maskR & tmat(:,l))-meandat)/stddat;
                volvecSafe(jj,k,l+1) = nanmean(hmat_PO(maskSafe & tmat(:,l))-meandat)/stddat;
                volvecRisky(jj,k,l+1) = nanmean(hmat_PO(maskRisky & tmat(:,l))-meandat)/stddat;
            end
        else
            volvecChoose(jj,k,1) = nanmean(hmat_reward(mask))/nf;
            volvecL(jj,k,1) = nanmean(hmat_PO(maskL))/nf;
            volvecR(jj,k,1) = nanmean(hmat_PO(maskR))/nf;
            %sem
            volvecChoose_sem(jj,k) = nanstd(hmat_reward(mask)/nf)/sqrt(sum(mask));
            volvecL_sem(jj,k) = nanstd(hmat_PO(maskL)/nf)/sqrt(sum(maskL));
            volvecR_sem(jj,k) = nanstd(hmat_PO(maskR)/nf)/sqrt(sum(maskR));

            for l = 1:nsamp
                volvecChoose(jj,k,l+1) = nanmean(hmat_reward(mask & tmat(:,l)))/nf;
                volvecL(jj,k,l+1) = nanmean(hmat_PO(maskL & tmat(:,l)))/nf;
                volvecR(jj,k,l+1) = nanmean(hmat_PO(maskR & tmat(:,l)))/nf;
            end
        end

    end

    %build the prob feature
    for k = 1:nbase_prob      

        mask = A{j}.chosenprob==probvec(k); %chosen value
        maskL = A{j}.left_prob==probvec(k);
        maskR = A{j}.right_prob==probvec(k);

        maskSafe =  ( A{j}.left_prob==probvec(k) & safeOptL) | ...
                    (A{j}.right_prob==probvec(k) & safeOptR);
        maskRisky = (A{j}.left_prob==probvec(k) & riskyOptL) | ...
                    (A{j}.right_prob==probvec(k) & riskyOptR);   

        if useZscore
            probvecChoose(jj,k,1) = nanmean(hmat_reward(mask)-meandat)/stddat;
            probvecL(jj,k,1) = nanmean(hmat_PO(maskL)-meandat)/stddat;
            probvecR(jj,k,1) = nanmean(hmat_PO(maskR)-meandat)/stddat;
            probvecSafe(jj,k,1) = nanmean(hmat_PO(maskSafe)-meandat)/stddat;
            probvecRisky(jj,k,1) = nanmean(hmat_PO(maskRisky)-meandat)/stddat;
            %sem
            probvecChoose_sem(jj,k) = nanstd((hmat_reward(mask)-meandat)/stddat)/sqrt(sum(mask));
            probvecL_sem(jj,k) = nanstd((hmat_PO(maskL)-meandat)/stddat)/sqrt(sum(maskL));
            probvecR_sem(jj,k) = nanstd((hmat_PO(maskR)-meandat)/stddat)/sqrt(sum(maskR));
            probvecSafe_sem(jj,k) = nanstd((hmat_PO(maskSafe)-meandat)/stddat)/sqrt(sum(maskSafe));
            probvecRisky_sem(jj,k) = nanstd((hmat_PO(maskRisky)-meandat)/stddat)/sqrt(sum(maskRisky));

            for l = 1:nsamp
                probvecChoose(jj,k,l+1) = nanmean(hmat_reward(mask & tmat(:,l))-meandat)/stddat;
                probvecL(jj,k,l+1) = nanmean(hmat_PO(maskL & tmat(:,l))-meandat)/stddat;
                probvecR(jj,k,l+1) = nanmean(hmat_PO(maskR & tmat(:,l))-meandat)/stddat;                
                probvecSafe(jj,k,l+1) = nanmean(hmat_PO(maskSafe & tmat(:,l))-meandat)/stddat;
                probvecRisky(jj,k,l+1) = nanmean(hmat_PO(maskRisky & tmat(:,l))-meandat)/stddat;
            end

        else
            probvecChoose(jj,k,1) = nanmean(hmat_reward(mask))/nf;
            probvecL(jj,k,1) = nanmean(hmat_PO(maskL))/nf;
            probvecR(jj,k,1) = nanmean(hmat_PO(maskR))/nf;
            %sem
            probvecChoose_sem(jj,k) = nanstd(hmat_reward(mask)/nf)/sqrt(sum(mask));
            probvecL_sem(jj,k) = nanstd(hmat_PO(maskL)/nf)/sqrt(sum(maskL));
            probvecR_sem(jj,k) = nanstd(hmat_PO(maskR)/nf)/sqrt(sum(maskR));

            for l = 1:nsamp
                probvecChoose(jj,k,l+1) = nanmean(hmat_reward(mask & tmat(:,l)))/nf;
                probvecL(jj,k,l+1) = nanmean(hmat_PO(maskL & tmat(:,l)))/nf;
                probvecR(jj,k,l+1) = nanmean(hmat_PO(maskR & tmat(:,l)))/nf;
            end
        end

    end
    
    %rewarded volume features
    rewvol = A{j}.chosenval.*A{j}.hits;    
    
    for k = 1:nbase_rew
        maskrewVol = rewvol==EV_bins(k);        
        if useZscore
            rewVolvec(jj,k,1) = nanmean(hmat_reward(maskrewVol)-meandat)/stddat;
            rewVolvec_sem(jj,k) = nanstd((hmat_reward(maskrewVol)-meandat)/stddat)/sqrt(sum(maskrewVol));           
            for l = 1:nsamp
                rewVolvec(jj,k,l+1) = nanmean(hmat_reward(maskrewVol & tmat(:,l))-meandat)/stddat;                
            end          
        else
             rewVolvec(jj,k,1) = nanmean(hmat_reward(maskrewVol))/nf;
             rewVolvec_sem(jj,k) = nanstd(hmat_reward(maskrewVol)/nf)/sqrt(sum(maskrewVol));
             
            for l = 1:nsamp
                rewVolvec(jj,k,l+1) = nanmean(hmat_reward(maskrewVol & tmat(:,l)))/nf;               
            end          
        end      
    end
    


    %build the EV feature vec, and rewarded value
    ev = A{j}.chosenprob.*A{j}.chosenval;
    evL = A{j}.left_prob.*A{j}.this_left_volume;
    evR = A{j}.right_prob.*A{j}.this_right_volume;
    

    ev_binned = discretize(ev,EV_bins);
    evL_binned = discretize(evL,EV_bins);
    evR_binned = discretize(evR,EV_bins);

    for k = 1:nbase_EV      

        mask = ev_binned==k; %chosen value
        maskL = evL_binned==k;
        maskR = evR_binned==k;
        

        maskSafe =  ( evL_binned==k & safeOptL) | ...
                    ( evR_binned==k & safeOptR);
        maskRisky = (evL_binned==k & riskyOptL) | ...
                    (evR_binned==k & riskyOptR);  


        if useZscore          
            EVvecChoose(jj,k,1) = nanmean(hmat_reward(mask)-meandat)/stddat;
            EVvecL(jj,k,1) = nanmean(hmat_PO(maskL)-meandat)/stddat;
            EVvecR(jj,k,1) = nanmean(hmat_PO(maskR)-meandat)/stddat;
            %use left-right hmats. they are what we want
            EVvecSafe(jj,k,1) = nanmean(hmat_PO(maskSafe)-meandat)/stddat;
            EVvecRisky(jj,k,1) = nanmean(hmat_PO(maskRisky)-meandat)/stddat;
            %sem          
            EVvecChoose_sem(jj,k) = nanstd((hmat_reward(mask)-meandat)/stddat)/sqrt(sum(mask));
            EVvecL_sem(jj,k) = nanstd((hmat_PO(maskL)-meandat)/stddat)/sqrt(sum(maskL));
            EVvecR_sem(jj,k) = nanstd((hmat_PO(maskR)-meandat)/stddat)/sqrt(sum(maskR));
            EVvecSafe_sem(jj,k) = nanstd((hmat_PO(maskSafe)-meandat)/stddat)/sqrt(sum(maskSafe));
            EVvecRisky_sem(jj,k) = nanstd((hmat_PO(maskRisky)-meandat)/stddat)/sqrt(sum(maskRisky));

            for l = 1:nsamp
                EVvecChoose(jj,k,l+1) = nanmean(hmat_reward(mask & tmat(:,l))-meandat)/stddat;
                EVvecL(jj,k,l+1) = nanmean(hmat_PO(maskL & tmat(:,l))-meandat)/stddat;
                EVvecR(jj,k,l+1) = nanmean(hmat_PO(maskR & tmat(:,l))-meandat)/stddat;
                EVvecSafe(jj,k,l+1) = nanmean(hmat_PO(maskSafe & tmat(:,l))-meandat)/stddat;
                EVvecRisky(jj,k,l+1) = nanmean(hmat_PO(maskRisky & tmat(:,l))-meandat)/stddat;
            end
        else   
            EVvecChoose(jj,k,1) = nanmean(hmat_reward(mask))/nf;
            EVvecL(jj,k,1) = nanmean(hmat_PO(maskL))/nf;
            EVvecR(jj,k,1) = nanmean(hmat_PO(maskR))/nf;
            %sem
            EVvecChoose_sem(jj,k) = nanstd(hmat_reward(mask)/nf)/sqrt(sum(mask));
            EVvecL_sem(jj,k) = nanstd(hmat_PO(maskL)/nf)/sqrt(sum(maskL));
            EVvecR_sem(jj,k) = nanstd(hmat_PO(maskR)/nf)/sqrt(sum(maskR));

            for l = 1:nsamp               
                EVvecChoose(jj,k,l+1) = nanmean(hmat_reward(mask & tmat(:,l)))/nf;
                EVvecL(jj,k,l+1) = nanmean(hmat_PO(maskL & tmat(:,l)))/nf;
                EVvecR(jj,k,l+1) = nanmean(hmat_PO(maskR & tmat(:,l)))/nf;
            end
        end

    end  


    %RPE adn delta EV, and rewarded volume
    ev = A{j}.chosenprob.*A{j}.chosenval;
    rewval = A{j}.chosenval.*A{j}.hits;
    rpe = rewval-ev;
    rpe_binned = discretize(rpe,RPE_bins);
    dev = evL-evR;
    dev_binned =  discretize(dev,DEV_bins);

    for k = 1:nbase_RPE

        mask = rpe_binned==k; 
        maskDEV = dev_binned==k; 


        if useZscore
            RPEvec(jj,k,1) = nanmean(hmat_reward(mask)-meandat)/stddat;
            RPEvec_sem(jj,k) = nanstd((hmat_reward(mask)-meandat)/stddat)/sqrt(sum(mask));

            DEVvec(jj,k,1) = nanmean(hmat_PO(maskDEV)-meandat)/stddat;
            DEVvec_sem(jj,k) = nanstd((hmat_PO(maskDEV)-meandat)/stddat)/sqrt(sum(maskDEV));

            for l = 1:nsamp
                RPEvec(jj,k,l+1) = nanmean(hmat_reward(mask & tmat(:,l))-meandat)/stddat;
                DEVvec(jj,k,l+1) = nanmean(hmat_PO(maskDEV & tmat(:,l))-meandat)/stddat;
            end

        else
            RPEvec(jj,k,1) = nanmean(hmat_reward(mask))/nf;
            RPEvec_sem(jj,k) = nanstd(hmat_reward(mask)/nf)/sqrt(sum(mask));

            DEVvec(jj,k,1) = nanmean(hmat_PO(maskDEV))/nf;
            DEVvec_sem(jj,k) = nanstd(hmat_PO(maskDEV)/nf)/sqrt(sum(maskDEV));

            for l = 1:nsamp
                RPEvec(jj,k,l+1) = nanmean(hmat_reward(mask & tmat(:,l)))/nf;
                DEVvec(jj,k,l+1) = nanmean(hmat_PO(maskDEV & tmat(:,l)))/nf;
            end
        end


    end

end

%add to the correct cells

%Rewarded value offer features
Rewfet_legend = {'volChoose','probChoose','EVChoose','RPE','rewVol'};
Rewfet = {volvecChoose,probvecChoose, EVvecChoose, RPEvec,rewVolvec};
Rewfet_sem = {volvecChoose_sem,probvecChoose_sem, EVvecChoose_sem, RPEvec_sem,rewVolvec_sem};

%presented value offer features
POfet_legend = {'volL','volR','probL','probR','EVL','EVR','DEV',...
    'volRisky','volSafe','probRisky','probSafe','EVRisky','EVSafe'};
POfet = {volvecL,volvecR, probvecL, probvecR, EVvecL, EVvecR, DEVvec,...
    volvecRisky, volvecSafe, probvecSafe, probvecRisky, EVvecRisky, EVvecSafe};
POfet_sem = {volvecL_sem,volvecR_sem, probvecL_sem, probvecR_sem,...
    EVvecL_sem, EVvecR_sem, DEVvec_sem,...
    volvecRisky_sem, volvecSafe_sem, probvecSafe_sem,...
    probvecRisky_sem, EVvecRisky_sem, EVvecSafe_sem};


save(strcat(ops.savedir,fname),'COfet','COfet_sem', 'COfet_legend',...
    'Rewfet','Rewfet_sem','Rewfet_legend',...
    'POfet','POfet_sem','POfet_legend',...
    'probvec','volvec','EV_bins','RPE_bins','DEV_bins',...
    'ops','sourceops');

%% create the conditional PSTHs of interest for comparisons

%TODO: include sem?

tmesh_start = -2:0.05:4;
tmesh_choice = -4:0.05:4;
nts = numel(tmesh_start);
ntc = numel(tmesh_choice);

%last ind=1 is regular, ind=2 is z-scored
psth = zeros(nn,nts,2); %fully marginalized
psth_choice = zeros(nn,ntc,2); %fully marginalized, alighned to choice
psth_leavecpoke = zeros(nn,ntc,2); %fully marginalized, align to leavecpoke
psth_prevWin = zeros(nn,nts,2);
psth_prevLoss = zeros(nn,nts,2);
psth_win = zeros(nn,ntc,2);
psth_loss = zeros(nn,ntc,2);
psth_left = zeros(nn,ntc,2);
psth_right = zeros(nn,ntc,2);
psth_risky = zeros(nn,ntc,2);
psth_safe = zeros(nn,ntc,2);

%expected value, volume, and probability
psth_evL = zeros(nn,ntc,nbase_EV,2);
psth_evR = zeros(nn,ntc,nbase_EV,2);
psth_volL = zeros(nn,ntc,nbase_vol,2);
psth_volR = zeros(nn,ntc,nbase_vol,2);
psth_probL = zeros(nn,ntc,nbase_prob,2);
psth_probR = zeros(nn,ntc,nbase_prob,2);
psth_probchosen = zeros(nn,ntc,nbase_prob,2);
%and for safe/risky
psth_evSafe = zeros(nn,ntc,nbase_EV,2);
psth_evRisky = zeros(nn,ntc,nbase_EV,2);
psth_volSafe = zeros(nn,ntc,nbase_vol,2);
psth_volRisky = zeros(nn,ntc,nbase_vol,2);
psth_probSafe = zeros(nn,ntc,nbase_prob,2);
psth_probRisky = zeros(nn,ntc,nbase_prob,2);

%rewarded volume, RPE
psth_rewVol = zeros(nn,ntc,nbase_rew,2);
psth_rpe = zeros(nn,ntc,nbase_RPE,2);
psth_dev = zeros(nn,ntc,nbase_RPE,2);

%the coarse reward volumes and expected volumes
psth_rewVol_coarse = zeros(nn,ntc,nbase_rew_coarse,2);
psth_rewVol_coarse_prevstart = zeros(nn,nts,nbase_rew_coarse,2); %same, but start aligned, and previous reward
psth_ev_chosen_coarse = zeros(nn,ntc,nbase_EV_coarse,2); %expected value of chosen option
psth_ev_chosen_prevstart_coarse = zeros(nn,nts,nbase_EV_coarse,2); %expected value of chosen option,start aligned

%the preference-parametrized responses. var1 is higher response than var2
psth_pref_prevWinLoss1 = zeros(nn,nts,2);
psth_pref_prevWinLoss2 = zeros(nn,nts,2);
psth_pref_winLoss1 = zeros(nn,ntc,2);
psth_pref_winLoss2 = zeros(nn,ntc,2);
psth_pref_leftRight1 = zeros(nn,ntc,2);
psth_pref_leftRight2 = zeros(nn,ntc,2);
psth_pref_riskySafe1 = zeros(nn,ntc,2);
psth_pref_riskySafe2 = zeros(nn,ntc,2);

disp('building cond. psths for later visualizations')

for jj = 1:nn
    disp(j) 
    j = Aj(jj); 
   
    %get max firing rate of cell to normalize
    %nf = nanmax(A{j}.hmat_start,[],'all');
    meandat = nanmean(nanmean(A{j}.hmat_start)); %unsure if I should use mean of mean, or mean of all?
    stddat = nanstd(nanmean(A{j}.hmat_start));
    
    meandat_choice = nanmean(nanmean(A{j}.hmat_choice));
    stddat_choice = nanstd(nanmean(A{j}.hmat_choice));
    meandat_leavecpoke = nanmean(nanmean(A{j}.hmat_leavecpoke));
    stddat_leavecpoke = nanstd(nanmean(A{j}.hmat_leavecpoke));
    
    dat1 = A{j}.hmat_start;
    dat2 = A{j}.hmat_choice;
    dat3 = A{j}.hmat_leavecpoke;
    
    %the full psth
    psth(jj,:,1) =  nanmean(dat1);
    psth(jj,:,2) = (nanmean(dat1)-meandat)/stddat;
    
    psth_choice(jj,:,1) =  nanmean(dat2);
    psth_choice(jj,:,2) = (nanmean(dat2)-meandat_choice)/stddat_choice;
    psth_leavecpoke(jj,:,1) =  nanmean(dat3);
    psth_leavecpoke(jj,:,2) = (nanmean(dat3)-meandat_leavecpoke)/stddat_leavecpoke;
    
    %the binary conditional PSTHs
    psth_win(jj,:,1) = nanmean(dat2(mask_cell{1}(j),:));
    psth_win(jj,:,2) = (nanmean(dat2(mask_cell{1}(j),:))-meandat)/stddat;
    
    psth_loss(jj,:,1) = nanmean(dat2(mask_cell{2}(j),:));
    psth_loss(jj,:,2) = (nanmean(dat2(mask_cell{2}(j),:))-meandat)/stddat;
    
    psth_prevWin(jj,:,1) = nanmean(dat1(mask_cell{3}(j),:));
    psth_prevWin(jj,:,2) = (nanmean(dat1(mask_cell{3}(j),:))-meandat)/stddat;
    
    psth_prevLoss(jj,:,1) = nanmean(dat1(mask_cell{4}(j),:));
    psth_prevLoss(jj,:,2) = (nanmean(dat1(mask_cell{4}(j),:))-meandat)/stddat;
    
    psth_left(jj,:,1) = nanmean(dat3(mask_cell{5}(j),:));
    psth_left(jj,:,2) = (nanmean(dat3(mask_cell{5}(j),:))-meandat)/stddat;
    
    psth_right(jj,:,1) = nanmean(dat3(mask_cell{6}(j),:));
    psth_right(jj,:,2) = (nanmean(dat3(mask_cell{6}(j),:))-meandat)/stddat;
    
    psth_risky(jj,:,1) = nanmean(dat3(mask_cell{7}(j),:));
    psth_risky(jj,:,2) = (nanmean(dat3(mask_cell{7}(j),:))-meandat)/stddat;
    
    psth_safe(jj,:,1) = nanmean(dat3(mask_cell{8}(j),:));
    psth_safe(jj,:,2) = (nanmean(dat3(mask_cell{8}(j),:))-meandat)/stddat;
    
    
    
    %the "preference" responses
    
    %win vs. loss
    tmesh = -4:0.05:4;
    wind_idx = discretize(wind_CO{1},tmesh); %win-loss indices
    c1 = nanmean(nanmean(dat2(mask_cell{1}(j),wind_idx(1):wind_idx(2)),2)); %win
    c2 = nanmean(nanmean(dat2(mask_cell{2}(j),wind_idx(1):wind_idx(2)),2)); %loss
    [~,pref_winloss] = min([c1,c2]);
    if pref_winloss==1 %win > loss
        psth_pref_winLoss1(jj,:,1) = nanmean(dat2(mask_cell{1}(j),:));
        psth_pref_winLoss1(jj,:,2) = (nanmean(dat2(mask_cell{1}(j),:))-meandat)/stddat;
        
        psth_pref_winLoss2(jj,:,1) = nanmean(dat2(mask_cell{2}(j),:));
        psth_pref_winLoss2(jj,:,2) = (nanmean(dat2(mask_cell{2}(j),:))-meandat)/stddat;
    else %loss > win
        psth_pref_winLoss1(jj,:,1) = nanmean(dat2(mask_cell{2}(j),:));
        psth_pref_winLoss1(jj,:,2) = (nanmean(dat2(mask_cell{2}(j),:))-meandat)/stddat;
        
        psth_pref_winLoss2(jj,:,1) = nanmean(dat2(mask_cell{1}(j),:));
        psth_pref_winLoss2(jj,:,2) = (nanmean(dat2(mask_cell{1}(j),:))-meandat)/stddat;
    end
    

    %prevWin vs. prevLoss
    %win vs. loss
    tmesh = -2:0.05:4;
    wind_idx = discretize(wind_CO{3},tmesh); %win-loss indices
    c1 = nanmean(nanmean(dat1(mask_cell{3}(j),wind_idx(1):wind_idx(2)),2)); %prevwin
    c2 = nanmean(nanmean(dat1(mask_cell{4}(j),wind_idx(1):wind_idx(2)),2)); %prevloss
    [~,pref_prevwinloss] = min([c1,c2]);
    if pref_prevwinloss ==1 %prevWin > prevLoss
        psth_pref_prevWinLoss1(jj,:,1) = nanmean(dat1(mask_cell{3}(j),:));
        psth_pref_prevWinLoss1(jj,:,2) = (nanmean(dat1(mask_cell{3}(j),:))-meandat)/stddat;
        
        psth_pref_prevWinLoss2(jj,:,1) = nanmean(dat1(mask_cell{4}(j),:));
        psth_pref_prevWinLoss2(jj,:,2) = (nanmean(dat1(mask_cell{4}(j),:))-meandat)/stddat;
        
    else
        psth_pref_prevWinLoss1(jj,:,1) = nanmean(dat1(mask_cell{4}(j),:));
        psth_pref_prevWinLoss1(jj,:,2) = (nanmean(dat1(mask_cell{4}(j),:))-meandat)/stddat;
        
        psth_pref_prevWinLoss2(jj,:,1) = nanmean(dat1(mask_cell{3}(j),:));
        psth_pref_prevWinLoss2(jj,:,2) = (nanmean(dat1(mask_cell{3}(j),:))-meandat)/stddat;
        
    end
    
    %left vs. right
    tmesh = -4:0.05:4;
    wind_idx = discretize(wind_CO{5},tmesh); %win-loss indices
    c1 = nanmean(nanmean(dat2(mask_cell{5}(j),wind_idx(1):wind_idx(2)),2)); %left
    c2 = nanmean(nanmean(dat2(mask_cell{6}(j),wind_idx(1):wind_idx(2)),2)); %righjt
    [~,pref_leftright] = min([c1,c2]);
    if pref_leftright==1 %left > right
        psth_pref_leftRight1(jj,:,1) = nanmean(dat3(mask_cell{5}(j),:));
        psth_pref_leftRight1(jj,:,2) = (nanmean(dat3(mask_cell{5}(j),:))-meandat)/stddat;
        
        psth_pref_leftRight2(jj,:,1) = nanmean(dat3(mask_cell{6}(j),:));
        psth_pref_leftRight2(jj,:,2) = (nanmean(dat3(mask_cell{6}(j),:))-meandat)/stddat;
    else %right > left
        psth_pref_leftRight1(jj,:,1) = nanmean(dat3(mask_cell{6}(j),:));
        psth_pref_leftRight1(jj,:,2) = (nanmean(dat3(mask_cell{6}(j),:))-meandat)/stddat;
        
        psth_pref_leftRight2(jj,:,1) = nanmean(dat3(mask_cell{5}(j),:));
        psth_pref_leftRight2(jj,:,2) = (nanmean(dat3(mask_cell{5}(j),:))-meandat)/stddat;
    end
    
    
    %risky vs. safe
    tmesh = -4:0.05:4;
    wind_idx = discretize(wind_CO{7},tmesh); %win-loss indices
    c1 = nanmean(nanmean(dat2(mask_cell{7}(j),wind_idx(1):wind_idx(2)),2)); %left
    c2 = nanmean(nanmean(dat2(mask_cell{8}(j),wind_idx(1):wind_idx(2)),2)); %righjt
    [~,pref_riskysafe] = min([c1,c2]);
    if pref_riskysafe==1 %left > right
        psth_pref_riskySafe1(jj,:,1) = nanmean(dat3(mask_cell{7}(j),:));
        psth_pref_riskySafe1(jj,:,2) = (nanmean(dat3(mask_cell{7}(j),:))-meandat)/stddat;
        
        psth_pref_riskySafe2(jj,:,1) = nanmean(dat3(mask_cell{8}(j),:));
        psth_pref_riskySafe2(jj,:,2) = (nanmean(dat3(mask_cell{8}(j),:))-meandat)/stddat;
    else %right > left
        psth_pref_riskySafe1(jj,:,1) = nanmean(dat3(mask_cell{8}(j),:));
        psth_pref_riskySafe1(jj,:,2) = (nanmean(dat3(mask_cell{8}(j),:))-meandat)/stddat;
        
        psth_pref_riskySafe2(jj,:,1) = nanmean(dat3(mask_cell{7}(j),:));
        psth_pref_riskySafe2(jj,:,2) = (nanmean(dat3(mask_cell{7}(j),:))-meandat)/stddat;
    end
    
    %---------
    
    %risky/safe masks
    safeOptL = A{j}.left_prob ==1;
    safeOptR = A{j}.right_prob==1;
    riskyOptL = A{j}.left_prob ~=1 ;
    riskyOptR = A{j}.right_prob ~=1;
    
    
    %presnted expected value on left and right 
    evL = A{j}.left_prob.*A{j}.this_left_volume;
    evR = A{j}.right_prob.*A{j}.this_right_volume;
    evL_binned = discretize(evL,EV_bins);
    evR_binned = discretize(evR,EV_bins);
    
    
    
    for k = 1:nbase_EV      
        maskL = evL_binned==k;
        maskR = evR_binned==k;
        
        maskRisky = (evL_binned==k & riskyOptL) | ...
                    (evR_binned==k & riskyOptR);
        maskSafe = (evL_binned==k & safeOptL) | ...
                    (evR_binned==k & safeOptR);
                
        maskchooseL = (evL_binned==k & mask_cell{5}(j));
        maskchooseR = (evR_binned==k & mask_cell{6}(j));
        maskchoose = maskchooseL | maskchooseR;
        
        psth_evL(jj,:,k,1) = nanmean(dat2(maskL,:));
        psth_evL(jj,:,k,2) = (nanmean(dat2(maskL,:))-meandat)/stddat;
        
        psth_evR(jj,:,k,1) = nanmean(dat2(maskR,:));
        psth_evR(jj,:,k,2) = (nanmean(dat2(maskR,:))-meandat)/stddat;
        
        psth_evRisky(jj,:,k,1) = nanmean(dat2(maskRisky,:));
        psth_evRisky(jj,:,k,2) = (nanmean(dat2(maskRisky,:))-meandat)/stddat;
        
        psth_evSafe(jj,:,k,1) = nanmean(dat2(maskSafe,:));
        psth_evSafe(jj,:,k,2) = (nanmean(dat2(maskSafe,:))-meandat)/stddat;
       
    end    
    
    
    %coarse expected value of chosen option   
    evL_binned_coarse = discretize(evL,EV_bins_coarse);
    evR_binned_coarse = discretize(evR,EV_bins_coarse);
    evL_binned_coarse_prev = discretize([nan;evL(1:end-1)],EV_bins_coarse);
    evR_binned_coarse_prev = discretize([nan;evR(1:end-1)],EV_bins_coarse);
    
    for k = 1:nbase_EV_coarse  
        maskchooseL = (evL_binned_coarse==k & mask_cell{5}(j));
        maskchooseR = (evR_binned_coarse==k & mask_cell{6}(j));
        maskchoose = maskchooseL | maskchooseR;
        
        maskchooseL_prev = (evL_binned_coarse_prev==k & mask_cell{5}(j));
        maskchooseR_prev = (evR_binned_coarse_prev==k & mask_cell{6}(j));
        maskchoose_prev = maskchooseL_prev | maskchooseR_prev;
        
        %ev of chosen option
        psth_ev_chosen_coarse(jj,:,k,1) = nanmean(dat2(maskchoose,:));
        psth_ev_chosen_coarse(jj,:,k,2) = (nanmean(dat2(maskchoose,:))-meandat)/stddat;
        
        psth_ev_chosen_prevstart_coarse(jj,:,k,1) = nanmean(dat1(maskchoose_prev,:));
        psth_ev_chosen_prevstart_coarse(jj,:,k,2) = (nanmean(dat1(maskchoose_prev,:))-meandat)/stddat;
        
        
    end    
    
    
    %presented volume
    for k = 1:nbase_vol
        maskL = A{j}.this_left_volume==volvec(k);
        maskR = A{j}.this_right_volume==volvec(k);
        
        maskRisky = (maskL & riskyOptL) | ...
                    (maskR & riskyOptR);
        maskSafe = (maskL & safeOptL) | ...
                    (maskR & safeOptR);
        
        psth_volL(jj,:,k,1) = nanmean(dat2(maskL,:));
        psth_volL(jj,:,k,2) = (nanmean(dat2(maskL,:))-meandat)/stddat;
        
        psth_volR(jj,:,k,1) = nanmean(dat2(maskR,:));
        psth_volR(jj,:,k,2) = (nanmean(dat2(maskR,:))-meandat)/stddat; 
        
        psth_volSafe(jj,:,k,1) = nanmean(dat2(maskSafe,:));
        psth_volSafe(jj,:,k,2) = (nanmean(dat2(maskSafe,:))-meandat)/stddat;
        
        psth_volRisky(jj,:,k,1) = nanmean(dat2(maskRisky,:));
        psth_volRisky(jj,:,k,2) = (nanmean(dat2(maskRisky,:))-meandat)/stddat; 
    end
    
    %presented probability
    for k = 1:nbase_prob
        maskL = A{j}.left_prob==probvec(k);
        maskR = A{j}.right_prob==probvec(k);
        
        maskRisky = (maskL & riskyOptL) | ...
                    (maskR & riskyOptR);
        maskSafe = (maskL & safeOptL) | ...
                    (maskR & safeOptR);
                
        maskchosen = (A{j}.right_prob==probvec(k) & A{j}.went_right==1) |...
                        (A{j}.left_prob==probvec(k) & A{j}.went_right==0);
                    
        psth_probL(jj,:,k,1) = nanmean(dat2(maskL,:));
        psth_probL(jj,:,k,2) = (nanmean(dat2(maskL,:))-meandat)/stddat;
        
        psth_probR(jj,:,k,1) = nanmean(dat2(maskR,:));
        psth_probR(jj,:,k,2) = (nanmean(dat2(maskR,:))-meandat)/stddat;  
        
        psth_probSafe(jj,:,k,1) = nanmean(dat2(maskSafe,:));
        psth_probSafe(jj,:,k,2) = (nanmean(dat2(maskSafe,:))-meandat)/stddat;
        
        psth_probRisky(jj,:,k,1) = nanmean(dat2(maskRisky,:));
        psth_probRisky(jj,:,k,2) = (nanmean(dat2(maskRisky,:))-meandat)/stddat;  
        
        psth_probchosen(jj,:,k,1) = nanmean(dat2(maskchosen,:));
        psth_probchosen(jj,:,k,2) = (nanmean(dat2(maskchosen,:))-meandat)/stddat;  
    
    end
    
    %rewarded volume
    rewVol = A{j}.chosenval.*A{j}.hits;   
    for k = 1:nbase_rew
        mask = rewVol==EV_bins(k);       
        psth_rewVol(jj,:,k,1) = nanmean(dat2(mask,:));
        psth_rewVol(jj,:,k,2) = (nanmean(dat2(mask,:))-meandat)/stddat;
        
    end
    
    %coarse rewarded volume
    %rewarded volume
    rewVol = A{j}.chosenval.*A{j}.hits;   
    rewVol_prev = [nan;rewVol(1:end-1)];
    
    for k = 1:nbase_rew_coarse
        mask = rewVol==EV_bins_coarse(k);   
        maskprev = rewVol_prev==EV_bins_coarse(k);  
        
        psth_rewVol_coarse(jj,:,k,1) = nanmean(dat2(mask,:));
        psth_rewVol_coarse(jj,:,k,2) = (nanmean(dat2(mask,:))-meandat)/stddat;
        
        psth_rewVol_coarse_prevstart(jj,:,k,1) = nanmean(dat1(maskprev,:));
        psth_rewVol_coarse_prevstart(jj,:,k,2) = (nanmean(dat1(maskprev,:))-meandat)/stddat;
        
    end
    
    %RPE of reward
    ev = A{j}.chosenprob.*A{j}.chosenval;
    rpe = (A{j}.chosenval.*A{j}.hits)-ev;
    evL = A{j}.left_prob.*A{j}.this_left_volume;
    evR = A{j}.right_prob.*A{j}.this_right_volume;
    dev = evL-evR;
    rpe_binned = discretize(rpe,RPE_bins);
    dev_binned = discretize(dev,RPE_bins);
    
    for k = 1:nbase_RPE       
        mask = rpe_binned==k;       
        mask_dev = dev_binned==k;
        
        psth_rpe(jj,:,k,1) = nanmean(dat2(mask,:));
        psth_rpe(jj,:,k,2) = (nanmean(dat2(mask,:))-meandat)/stddat;
        
        psth_dev(jj,:,k,1) = nanmean(dat2(mask_dev,:));
        psth_dev(jj,:,k,2) = (nanmean(dat2(mask_dev,:))-meandat)/stddat;
    end
    
        
end


save(strcat(ops.savedir,'psth_conditional_all.mat'),'psth*','tmesh*',...
    'volvec','probvec','EV_bins*','RPE_bins','ops');


%% build the psth features

%the psth can be used as a feature. make a feature type object 
%include choice and leavecpoke aligned
psth_fet = zeros(nn,nts,nsamp+1);
psth_fet_sem = zeros(nn,nts);

psth_fet_choice = zeros(nn,ntc,nsamp+1);
psth_fet_choice_sem = zeros(nn,ntc);

psth_fet_leavecpoke = zeros(nn,ntc,nsamp+1);
psth_fet_leavecpoke_sem = zeros(nn,ntc);

tmesh = -2:0.05:4;
tmesh_choice = -4:0.05:4;
tmesh_leavecpoke = -4:0.05:4;

disp('building psth features')

for k = 1:3
    %decide which data to use
    switch k
        case 1
            datfun = @(j) A{j}.hmat_start; %function handle to get heatmat
            nt_k = nts; %time points
        case 2
            datfun = @(j) A{j}.hmat_choice;
            nt_k = ntc;
        case 3
            datfun = @(j) A{j}.hmat_leavecpoke;
            nt_k = ntc;
    end
        
    %the main variables. will be assinged at end of k loop
    psth_fet_k = zeros(nn,nt_k,nsamp+1);
    psth_fet_k_sem = zeros(nn,nt_k);

    for jj = 1:nn
        disp(jj)
        j = Aj(jj);
        
        dat = datfun(j);

        %get max firing rate of cell to normalize
            if normfet
                nf = nanmax(dat,[],'all');
            else
                nf = 1;
            end
            meandat = nanmean(nanmean(dat)); %unsure if I should use mean of mean, or mean of all?
            stddat = nanstd(nanmean(dat));

            %bootstrap masks
            ns = numel(A{j}.nspikes); %numbers total samples
            tmat = false(ns,nsamp);
            for m = 1:nsamp
                randoff = randsample(ns,floor(bootstrap_proportion*ns),bootstrap_replacemet);
                tmat(randoff,m) = true;
            end

            if useZscore
                psth_fet_k(jj,:,1) = nanmean(dat-meandat)/stddat;
                psth_fet_k_sem(jj,:) = nanstd((dat-meandat)/stddat)/sqrt(size(dat,2));

                for l = 1:nsamp
                    psth_fet_k(jj,:,l+1) = nanmean(dat(tmat(:,l),:)-meandat)/stddat;
                end
            else           
                psth_fet_k(jj,:,1) = nanmean(dat(mask))/nf;
                psth_fet_k_sem(jj,:) = nanstd(dat(mask)/nf)/sqrt(numel(dat(mask)));

                for l = 1:nsamp
                     psth_fet_k(jj,:,l+1) = nanmean(dat(tmat(:,l),:))/nf;
                end
            end
     
            
    end
    
    %assign psth_fet_k to correct output
    switch k
        case 1
            psth_fet = psth_fet_k;
            psth_fet_sem = psth_fet_k_sem;
        case 2
            psth_fet_choice = psth_fet_k;
            psth_fet_choice_sem = psth_fet_k_sem;
        case 3
            psth_fet_leavecpoke = psth_fet_k;
            psth_fet_leavecpoke_sem = psth_fet_k_sem;
    end
        
end

%append the source features to include the marginal PSTH as a feature type,
%with bootstraps
save(strcat(ops.savedir,fname),'psth_fet*','tmesh*','ops','-append');





