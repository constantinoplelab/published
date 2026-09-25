% parse up the beliefs into a 2d discrete distribution, but from raw matlab
% file

function [output] = decode_parsedata_rawdat(fname,varargin)
% parse up the beliefs into a 2d discrete distribution, but from raw matlab
% file. a general form of dSTR_decode_parsedata that uses Maggie's curated
% dSTR data. does beleifs only for a single session at a time, returns
% heatmats for all data
% INPUTS:
%   fname: filename of matlab file with S and SU struct
%   varargin: if present, (str) for which epoch to extract from heatmat
% OUTPUTS:
%   output: contains most likely block, 2nd mlb, x vector of timing, 
%           heatmat, etc.


%% general appraoch that should work for all cases

% for each neuron, get the following
% session_id
% neurons' average firing over time (2s after reward)
% p_mix and p_high for that trial
% which bin it is in

%% load code and neural/behavioral data
codepath = '~/projects/constantinoplelab/Analysis/';
% add code paths
addpath(genpath(codepath))
addpath(genpath(codepath + "david"))
addpath(genpath(codepath + "maggie"))

b = load(fname);

%% check if epoch was defined
if numel(varargin) == 1
    epoch = varargin{1};
else
    epoch = 'reward';
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

%% data used for regression


S = b.S;
SU = b.SU;

% this struct is almost the correct format for an A struct, but needs
% a little renaming to run. 
S.ntrials = [numel(S.Block)]; % only 1 session in this struct
S.reward = S.RewardAmount;
S.prob_catch = S.ProbCatch;

%neural data
%hmat_j = b.all_SU_DS{j}.hmat.Rew;
%xvec_j = b.all_SU_DS{j}.xvec.Rew;
% will need to downsample for 200ms bins
ncell = numel(SU);
hmats = cell(ncell,1);

%regen heatmats at 200ms precision, from -4 to 4s. to standardize hmat
bins = 0.2; %time bin size, in s
wndw = [-4,4];
behEvents = [];
if isfield(S,'behEvents')
    behEvents = S.behEvents; 
else
    behEvents = struct();
    behEvents.All = 1:length(S.NoseInCenter); 
    behEvents.Rewarded = find(S.hits); 
    behEvents.Optout = find(S.optout);
end
SU = makeHeatmat(SU, S, behEvents, wndw, bins, varargin);

for j = 1:ncell
% TODO: choose cell
    switch epoch
        case 'reward'
            hmats{j} = SU{j}.hmat.Rew;
            xvec = SU{1}.xvec.Rew;
        case 'coff'
            hmats{j} = SU{j}.hmat.COFF;
            xvec = SU{1}.xvec.COFF;
        case 'son'
            hmats{j} = SU{j}.hmat.SON;
            xvec = SU.xvec.SON;
        case 'con'
            hmats{j} = SU{j}.hmat.CON;
            xvec = SU{1}.xvec.CON;
    end
end

% Andrew's code:
% Analysis/BehavioralModel/anlaysis/BayesianModel/GenerateSynthData_Bayes.m
[~, ~, ~, Belief, ~, ~] = GenerateSynthData_Bayes(params, S,'logn', true,noise);

% get block and most likely block
[~,idx_max] = maxk(Belief,2,1);
idx_max(:,isnan(Belief(1,:))) = nan;
offer = S.reward;

mostlikely_block = idx_max(1,:);
secondlikely_block = idx_max(2,:);


%% todo: restructure
output = struct();
output.X = hmats;
output.xvec = xvec;
output.mostlikely_block = mostlikely_block;
output.secondlikely_block = secondlikely_block;
output.offer = offer;
output.rewarded_trials = S.hits==1;










