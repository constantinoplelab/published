% parse up the beliefs into a 2d discrete distribution

function [output] = dSTR_decode_parsedata(usetorch, varargin)


%% general appraoch that should work for all cases

% for each neuron, get the following
% session_id
% neurons' average firing over time (2s after reward)
% p_mix and p_high for that trial
% which bin it is in

%% load code data
codepath = '~/projects/constantinoplelab/Analysis/';

% add code paths
addpath(genpath(codepath))
addpath(genpath(codepath + "david"))
addpath(genpath(codepath + "maggie"))

% DS projection neuron data location
% change this: 
maggie = 1;
if usetorch
    datadir = "/scratch/dh148/dynamics/data/maggie/";
elseif maggie==0
    datadir = "/Users/dhocker/projects/dynamics/data/maggie/";
elseif maggie ==1
    datadir = '\\constantinoplelab.cns.nyu.edu\server\david\maggie\data\';
end

fname = datadir + "DS_projection_neurons_non-Stimulated.mat";
b = load(fname);

%% check if epoch was defined
if numel(varargin) == 1
    epoch = varargin{1};
    do_meansub = false;
elseif numel(varargin) == 2
    epoch = varargin{1};
    balance_blocks = varargin{2}; % will you balance blocks before marginalizing?
    if strcmp(balance_blocks,'none')
        do_meansub = false;
    else
        disp('doing mean subtraction of firing rate, conditioned on offer')
        do_meansub = true;
    end

else
    epoch = 'reward';
    do_meansub = false;
end

%% for every cell (only 60 now, and 23 sessions, so not prohibitive), do bayes calc

% bayes params for all sessions
kappa_mi = 1; %.5;  
kappa_hi = 1.5; %1.2; %0.7; 
kappa_lo = 0.5;%0.8; %.3; 
D = 0.5; 
lambda =1;
noise = 80;
params = [kappa_mi, kappa_hi, kappa_lo, D, lambda];

% data used for regression
ncell = size(b.all_index_DS,1);
sess_ids = zeros(ncell,1); % session id for this neuron
X_percell = cell(ncell,1); % neuron responses on each trial
mostlikely_block = cell(ncell,1); % the most likely block from each prob estimate
secondlikely_block = cell(ncell,1); % second most likely block
offer = cell(ncell,1); % offer associated with each belief
hits = cell(ncell,1);

%a concatenated version of targets and trials for better k-fold parsing
alltrials = []; %unique trial id across all sessions
all_mostlikey_block = []; % categories for all trials, concatenated
all_secondlikely_block = []; % second most likely 
all_hits = []; %same, for rewarded trial tracking
meanrates = [];

for j = 1:ncell
    % which session id for this cell?
    S_idx = b.all_index_DS{j,4}; 
    S_j = b.all_S_DS{S_idx};
    

    % this struct is almost the correct format for an A struct, but needs
    % a little renaming to run. 
    S_j.ntrials = [numel(S_j.Block)]; % only 1 session in this struct
    S_j.reward = S_j.RewardAmount;
    S_j.prob_catch = S_j.ProbCatch;

    %mean firing rate 
    meanrate_j = session_averaged_FR(b.all_SU_DLS_200{j}, 10);
    meanrates = [meanrates, meanrate_j];

    %neural data
    %hmat_j = b.all_SU_DS{j}.hmat.Rew;
    %xvec_j = b.all_SU_DS{j}.xvec.Rew;
    % 200ms bins instead
    switch epoch
        case 'reward'
            hmat_j = b.all_SU_DLS_200{j}.hmat.Rew;
            xvec_j = b.all_SU_DLS_200{j}.xvec.Rew;
        case 'coff'
            hmat_j = b.all_SU_DLS_200{j}.hmat.COFF;
            xvec_j = b.all_SU_DLS_200{j}.xvec.COFF;
        case 'son'
            hmat_j = b.all_SU_DLS_200{j}.hmat.SON;
            xvec_j = b.all_SU_DLS_200{j}.xvec.SON;
        case 'con'
            hmat_j = b.all_SU_DLS_200{j}.hmat.CON;
            xvec_j = b.all_SU_DLS_200{j}.xvec.CON;
    end

    if do_meansub
        switch balance_blocks
            case 'theoretical'

                disp('theoretical balancing of blocks:')
                fr_byvol_byblock = zeros(5,3);
                fr_byvol_balanced = zeros(5,1);
                vols = [5,10,20,40,80];
                for k = 1:5
                    for m = 1:3
                        mask = S_j.reward==vols(k) & S_j.Block == m;
                        fr_byvol_byblock(k,m) = mean(mean(hmat_j(mask,:),'omitnan'),'omitnan');            
                    end
    
                    % do empirical balancing
                    if k == 1 || k == 2
                        fr_byvol_balanced(k) = 0.5*fr_byvol_byblock(k,1) + 0.5*fr_byvol_byblock(k,3);
                    elseif k == 4 || k == 5
                        fr_byvol_balanced(k) = 0.5*fr_byvol_byblock(k,1) + 0.5*fr_byvol_byblock(k,2);
                    else
                        fr_byvol_balanced(k) = 1/3*fr_byvol_byblock(k,1) + 1/3*fr_byvol_byblock(k,2) + 1/3*fr_byvol_byblock(k,3);
                    end
    
    
                    hmat_j(mask,:) = hmat_j(mask,:)-fr_byvol_balanced(k);
                end

                disp('comparison of anlaytical balancing vs. raw marginalization')
                disp(fr_byvol_byblock)
                disp('...')
                disp(mean(fr_byvol_byblock,2,'omitnan'))
                disp('--')
                disp(fr_byvol_balanced)

            case 'psth'
                ntx = size(hmat_j,2);
                %disp('subtracting off psth by vol')
                % subtract the psth by vol
                fr_byvol_byblock = zeros(5,3,ntx);
                fr_byvol_balanced = zeros(5,1,ntx);
                vols = [5,10,20,40,80];
                for k = 1:5
                    for m = 1:3
                        mask = S_j.reward==vols(k) & S_j.Block == m;
                        fr_byvol_byblock(k,m,:) = mean(hmat_j(mask,:),'omitnan');            
                    end

                    % do empirical balancing
                    if k == 1 || k == 2
                        fr_byvol_balanced(k,:) = 0.5*fr_byvol_byblock(k,1,:) + 0.5*fr_byvol_byblock(k,3,:);
                    elseif k == 4 || k == 5
                        fr_byvol_balanced(k,:) = 0.5*fr_byvol_byblock(k,1,:) + 0.5*fr_byvol_byblock(k,2,:);
                    else
                        fr_byvol_balanced(k,:) = 1/3*fr_byvol_byblock(k,1,:) + 1/3*fr_byvol_byblock(k,2,:) + 1/3*fr_byvol_byblock(k,3,:);
                    end
    
    
                    hmat_j(mask,:) = hmat_j(mask,:)-fr_byvol_balanced(k,:);
                end

            case 'psth2' % i caught a bug in the code, so the fastest, safest way to rerun it was to make a new case
                ntx = size(hmat_j,2);
                %disp('subtracting off psth by vol')
                % subtract the psth by vol
                fr_byvol_byblock = zeros(5,3,ntx);
                fr_byvol_balanced = zeros(5,1,ntx);
                vols = [5,10,20,40,80];
                for k = 1:5
                    for m = 1:3
                        %mask = S_j.reward==vols(k) & S_j.Block == m;
                        mask = S_j.reward==vols(k) & S_j.Block == m & S_j.hits == 1;
                        %mask = S_j.reward==vols(k) & S_j.Block == m & S_j.vios == 0;
                        fr_byvol_byblock(k,m,:) = mean(hmat_j(mask,:),'omitnan');            
                    end

                    %mask_k = S_j.reward==vols(k);
                    mask_k = S_j.reward==vols(k) & S_j.hits == 1;
                    %mask_k = S_j.reward==vols(k) & S_j.vios == 0;
    
                    % do empirical balancing
                    if k == 1 || k == 2
                        fr_byvol_balanced(k,:) = 0.5*fr_byvol_byblock(k,1,:) + 0.5*fr_byvol_byblock(k,3,:);
                    elseif k == 4 || k == 5
                        fr_byvol_balanced(k,:) = 0.5*fr_byvol_byblock(k,1,:) + 0.5*fr_byvol_byblock(k,2,:);
                    else
                        fr_byvol_balanced(k,:) = 1/3*fr_byvol_byblock(k,1,:) + 1/3*fr_byvol_byblock(k,2,:) + 1/3*fr_byvol_byblock(k,3,:);
                    end
    
    
                    hmat_j(mask_k,:) = hmat_j(mask_k,:)-fr_byvol_balanced(k,:);
                    %hmat_j(mask,:) = hmat_j(mask,:)-fr_byvol_balanced(k,:);
                end


            otherwise

                disp('reward subtract, no block balancing')
                fr_byvol = zeros(5,1);
                vols = [5,10,20,40,80];
                for k = 1:5
                    mask = S_j.reward==vols(k);
                    fr_byvol(k) = mean(mean(hmat_j(mask,:),'omitnan'),'omitnan')
                    hmat_j(mask,:) = hmat_j(mask,:)-fr_byvol(k);
                end
        end
    end



    [~, wait_time, ~, Belief, ~, ~] = GenerateSynthData_Bayes(params, S_j,'logn', true,noise);

    % the long concatenated list of trial ids. encode session in there too.
    newtrials = 1:numel(wait_time);
    newtrials = newtrials + 10000*S_idx; % session identifier added
    alltrials = [alltrials, newtrials];

    % get block and most likely block
    [~,idx_max] = maxk(Belief,2,1);
    idx_max(:,isnan(Belief(1,:))) = nan;
    offer{j} = S_j.reward;
    hits{j} = S_j.hits';

    mostlikely_block{j} = idx_max(1,:);
    secondlikely_block{j} = idx_max(2,:);
    all_mostlikey_block = [all_mostlikey_block, mostlikely_block{j}];
    all_secondlikely_block = [all_secondlikely_block, secondlikely_block{j}];
    all_hits = [all_hits, S_j.hits'];

    sess_ids(j) = S_idx;
    X_percell{j} = hmat_j;

    

end

output = struct();
output.X = X_percell;
output.sess_ids = sess_ids;
output.xvec = xvec_j;
output.mostlikely_block = mostlikely_block;
output.secondlikely_block = secondlikely_block;
output.offer = offer;
output.meanrates = meanrates;
output.hits = hits;


% the concatenated ones
output.alltrials = alltrials;
output.all_mostlikey_block = all_mostlikey_block;
output.all_secondlikely_block = all_secondlikely_block;
output.all_hits = all_hits;








