% parse up the beliefs into a 2d discrete distribution

function [output] = dSTR_decode_parsedata(varargin)


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
datadir = "Z:\david\maggie\";
fname = datadir + "DS_projection_neurons_non-Stimulated.mat";
b = load(fname);

%% check if epoch was defined
if numel(varargin) == 1
    epoch = varargin{1}
else
    epoch = 'reward'
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

%a concatenated version of targets and trials for better k-fold parsing
alltrials = []; %unique trial id across all sessions
all_mostlikey_block = []; % categories for all trials, concatenated
all_secondlikely_block = []; % second most likely 
all_hits = []; %same, for rewarded trial tracking

for j = 1:ncell
    % which session id for this cell?
    S_idx = b.all_index_DS{j,4}; 
    S_j = b.all_S_DS{S_idx};

    % this struct is almost the correct format for an A struct, but needs
    % a little renaming to run. 
    S_j.ntrials = [numel(S_j.Block)]; % only 1 session in this struct
    S_j.reward = S_j.RewardAmount;
    S_j.prob_catch = S_j.ProbCatch;

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


    [~, wait_time, ~, Belief, ~, ~] = GenerateSynthData_Bayes(params, S_j,'logn', true,noise);

    % the long concatenated list of trial ids. encode session in there too.
    newtrials = 1:numel(wait_time);
    newtrials = newtrials + 10000*S_idx; % session identifier added
    alltrials = [alltrials, newtrials];

    % get block and most likely block
    [~,idx_max] = maxk(Belief,2,1);
    idx_max(:,isnan(Belief(1,:))) = nan;
    offer{j} = S_j.reward;

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


% the concatenated ones
output.alltrials = alltrials;
output.all_mostlikey_block = all_mostlikey_block;
output.all_secondlikely_block = all_secondlikely_block;
output.all_hits = all_hits;








