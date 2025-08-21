function ProcessData_Figure3(datadir, savedir, codedir)
%ProcessData_Figure3 - Process raw data saved under datadir such that it can be plotted by PlotFigure3.
% INPUTS:
%   datadir - Local directory where 'Golden_Nature/RawData_Figure3.mat' from Zenodo was saved
%   codedir - Local directory where code (e.g., published/EstrousRPE/) was saved
%   savedir - Local directory where you would like the outputs to be saved

%% Set up paths
s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codedir, 'IgnoreCase', ispc);
if ~onPath % only add code dir to path if it isn't already
    addpath(genpath(codedir))
end

%% Load raw data
load([datadir, 'RawData_Figure3'], 'ratTrial',...
    'BestFits', 'Pstruct', 'Bstruct', 'times_resampled', 'PstructCPIn',...
    'Proestrus_last', 'Estrus_last', 'Metestrus_last', 'Diestrus_last',...
    'pro_ratList', 'est_ratList', 'met_ratList', 'di_ratList')

%% Process data
%Set general variables
lat_thresh = 99;
Stages = {'Proestrus','Diestrus'}; %estrous stages to include
nback = 7;
NAcc_ratlist = {'G008','G016','G021','G022','G024',...
    'G036','G037','G051','G127'};

%--------------------------------------------------------------------------
%3a. Reinforcement learning model that computes initiation times from kappa
%(expected value). Example predicted initiation time and kappa over blocks.
%--------------------------------------------------------------------------
[~, ratTrial.reward] = convertreward(ratTrial.reward); %using standard reward amounts
%set first ITI of each session to NaN
ctr = 1;
for jk = 1:length(ratTrial.ntrials)
    ratTrial.ITI(ctr) = nan;
    ctr = ctr+ratTrial.ntrials(jk);
end
%remove outliers
ratTrial.ITI(ratTrial.ITI>prctile(ratTrial.ITI,lat_thresh)) = NaN;
alpha = 0.6;
opto_arg = 0;
gain_arg = 0;
fig = 3;
[LatOpt, ~, Kappa, ~] =...
    GenerateLatencyData_VanillaAlpha_cg(alpha, ratTrial, opto_arg,...
    gain_arg, fig);
%remove violation trials
Kappa_trials = Kappa(ratTrial.vios~=1);
Lat_trials = LatOpt(ratTrial.vios~=1);

%--------------------------------------------------------------------------
%3c. Regression of previous trial reward volumes on predicted initiation
%times
%--------------------------------------------------------------------------
ratModel = ratTrial;
ratModel.ITI = LatOpt;
[betas_model,~,~] =...
    regress_latency_vs_rew(ratModel, nback, 0, 0, 0);

%--------------------------------------------------------------------------
%3d. Regression coefficient by time window between model-predicted RPE and 
% trial-by-trial DA z-score (last 20 trials) across time bins around every
% event. Use 100 ms bins. Pooling  across sessions for each rat, average 
% over rats. Includes 0.5 s baseline before CPIn.
%--------------------------------------------------------------------------
%Fit model to pooled sessions across rats and use that RPE
Alignments = {'CPOn', 'CPIn', 'SideOn', 'SideOff', 'Reward', 'OptOut'};
EncodingWindowAnalysis = [];
for rat = 1:length(NAcc_ratlist)
    EncodingWindows = DARegressionAnalysis_CG(NAcc_ratlist(rat),...
        times_resampled, 0.1, Alignments, BestFits, Bstruct, Pstruct);
    EncodingWindowAnalysis.(NAcc_ratlist{rat}) = EncodingWindows;
end

%--------------------------------------------------------------------------
%3e. DA AUC at offer cue as a function of model-predicted (using alpha fit to the
% 10 trials) RPE
%--------------------------------------------------------------------------
numbins = 6;
window = 0.5; %AUC window
event = 'CPIn'; %align to reward offer cue
[pro_DA_binned, ~, ~, di_DA_binned,...
    RPEbins] = DA_by_RPE_estrous(NAcc_ratlist,...
    PstructCPIn, Bstruct, Proestrus_last, Estrus_last, Metestrus_last,...
    Diestrus_last, pro_ratList, est_ratList, met_ratList, di_ratList,...
    numbins, window, event);

%--------------------------------------------------------------------------
%3f-j, l-m. Predicted block sensitivity and regression with added RPE to larger rewards
%--------------------------------------------------------------------------
twin = 30;
smoothfactor = 15;
LatOpts = cell(2, 1);
RPEs = cell(2, 1);
RPEs_alphachange = cell(2, 1);
delta = cell(2, 1);
isblockchange = cell(2, 1);
ltom = cell(2, 1);
htom = cell(2, 1);
mtol = cell(2, 1);
mtoh = cell(2, 1);
betas_model_3j = cell(2, 1);
betas_RPE = cell(2, 1);
betas_RPE_alphachange = cell(2, 1);
for jj=1:2
    if jj==1 %proestrus model
        %don't change alpha, increase RPE gain factor
        params = [0.6 1.75]; %alpha, RPE gain factor
        gain_arg = 1;
        %increase alpha, don't change RPE gain factor
        params_alphachange = [0.8 1]; %alpha, RPE gain factor
        gain_arg_alphachange = 0;
    elseif jj==2 %diestrus model
        %don't change alpha, increase RPE gain factor
        params = [0.6 1]; %alpha, RPE gain factor
        gain_arg=1;       
        %increase alpha, don't change RPE gain factor
        params_alphachange = [0.6 1]; %alpha, RPE gain factor
        gain_arg_alphachange = 0;
    end
    %manipulating RPE gain
    [LatOpts{jj}, RPEs{jj}] =...
        GenerateLatencyData_VanillaAlpha_cg(params, ratTrial,... 
            0, gain_arg, 3);
    %manipulating alpha
    [~, RPEs_alphachange{jj}, ~, ~] =...
        GenerateLatencyData_VanillaAlpha_cg(params_alphachange, ratTrial,... 
            0, gain_arg_alphachange, 3);
    Block = ratTrial.block;
    RPE_alphachange = RPEs_alphachange{jj};
    iti_block_avg = NaN(1,2);
    L = LatOpts{jj};
    for bl = 2:3
        iti_block_avg(bl) = mean(L(Block==bl), 'omitnan');
    end
    delta{jj} = iti_block_avg(3) - iti_block_avg(2);

    %3g RPE over blocks of trials to show gain
    blockdifference = diff(Block);
    isblockchange{jj} = find(blockdifference~=0)+1;  

    % 3i latency dynamics
    ratModel = ratTrial;
    ratModel.ITI = LatOpts{jj};    
    [ltom{jj}, htom{jj}, mtol{jj}, mtoh{jj}] = block_dynamics_latency(ratModel,...
        twin, smoothfactor, 0, 1, 0, 1); %detrend and use causual smoothing

    %3j initiation time regression coefficients
    [betas_model_3j{jj},~,~] =...
        regress_latency_vs_rew(ratModel, nback, 0, 0, 0);

    %3l RPEs regression coefficients, RPE gain manipulation
    ratModelwRPE = ratTrial;
    ratModelwRPE.RPE = RPEs{jj};
    betas_RPE{jj} = regress_RPE_vs_rew(ratModelwRPE, nback);

    %3m RPEs regression coefficients, alpha manipulation
    ratModelwRPE.RPE = RPE_alphachange;
    betas_RPE_alphachange{jj} = regress_RPE_vs_rew(ratModelwRPE, nback);
    
end

%--------------------------------------------------------------------------
%3f Plot RPE as a function of gain applied to RPEs
%--------------------------------------------------------------------------
[ymean_RPEgain_binned, x_RPE_binned] = bin_y_for_x(RPEs{2}, RPEs{1}, 1);

%--------------------------------------------------------------------------
%3k. Regression coefficients of DA as a function of current and previous 
% rewards in mixed blocks
%--------------------------------------------------------------------------
auc1 = 0; %start of AUC window
auc2 = 0.5; %end of AUC window
[BetasPro, BetasDi] = DARegressRewardHistory_estrous(NAcc_ratlist,...
    PstructCPIn, Bstruct, auc1, auc2, nback, Stages);

%% Save processed data
save([savedir 'ProcessData_Figure3'],...
    'ratTrial', 'LatOpt', 'Kappa_trials', 'Lat_trials', 'betas_model',...
    'Alignments', 'NAcc_ratlist', 'EncodingWindowAnalysis',...
    'pro_DA_binned', 'di_DA_binned', 'RPEbins', 'LatOpts', 'RPEs',...
    'RPEs_alphachange', 'delta', 'isblockchange', 'ltom', 'htom',...
    'mtol', 'mtoh', 'betas_model_3j', 'betas_RPE',...
    'betas_RPE_alphachange', 'ymean_RPEgain_binned', 'x_RPE_binned',...
    'BetasPro', 'BetasDi');

end