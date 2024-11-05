function ProcessData_Figure1(datadir, codedir, savedir) 
%ProcessData_Figure1 - Process raw data saved under datadir such that it can be plotted by PlotFigure1.
% INPUTS:
%   datadir - Local directory where 'Golden_Nature/RawData_Figure1.mat' from Zenodo was saved
%   codedir - Local directory where code (e.g., published/EstrousRPE/) was saved
%   savedir - Local directory where you would like the outputs to be saved

%% Set up paths
s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codedir, 'IgnoreCase', ispc);
if ~onPath % only add code dir to path if it already isn't
    addpath(genpath(codedir))
end

%% Load raw data
%Rat lists
load([datadir, 'RawData_Figure1'], 'f_ratlist', 'm_ratlist',...
    'RatBehaviorData', 'SerumTable', 'SerumBehData');

%% Process data
%Initiation times by block
% Get data over estrous
ITIbyBlock_est = ITIbyBlock_estrous(f_ratlist, RatBehaviorData); 
ITIbyBlock.Proestrus = ITIbyBlock_est.Proestrus;
ITIbyBlock.Estrus = ITIbyBlock_est.Estrus;
ITIbyBlock.Metestrus = ITIbyBlock_est.Metestrus;
ITIbyBlock.Diestrus = ITIbyBlock_est.Diestrus;
beh_sens_staged = ((cell2mat(ITIbyBlock.Proestrus.low_ITI_det) + cell2mat(ITIbyBlock.Diestrus.low_ITI_det))./2) >...
    ((cell2mat(ITIbyBlock.Proestrus.high_ITI_det) + cell2mat(ITIbyBlock.Diestrus.high_ITI_det))./2); %require that low > high from this point on
new_f_ratlist = f_ratlist(beh_sens_staged);
newfrats = length(new_f_ratlist);

% Get data from females
ITIbyBlock_female = ITIbyBlock_females(f_ratlist, RatBehaviorData);
ITIbyBlock.Female = ITIbyBlock_female.Female;

% Get data from males
ITIbyBlock_male = ITIbyBlock_males(m_ratlist, RatBehaviorData); %not screening, excluding muscimol rats
ITIbyBlock.Male = ITIbyBlock_male.Male;
beh_sens_males = cell2mat(ITIbyBlock.Male.low_ITI_det) >...
    cell2mat(ITIbyBlock.Male.high_ITI_det); %require that low > high from this point on
new_m_ratlist = m_ratlist(beh_sens_males);
newmrats = length(new_m_ratlist);

%initiation times at block transitions
smoothfactor = 15;
twin = 30;
Latencydynamics = block_latencydynamics_macro_stages(twin,...
    smoothfactor,98,f_ratlist,RatBehaviorData);

% initiation times and estradiol levels from sessions with serum draws
detrend_arg = 1;
[delta_rats, E2_rats] =...
    Estradiol_vs_DeltaInit_bystage(SerumTable, SerumBehData, detrend_arg); 
%sessions were not screened for wait time by reward volume curves because
%only a subset were sampled

%% Save processed data
save([savedir 'ProcessData_Figure1'], 'f_ratlist', 'm_ratlist',...
    'newfrats', 'beh_sens_staged', 'newmrats', 'beh_sens_males',...
    'ITIbyBlock', 'Latencydynamics', 'delta_rats', 'E2_rats')

end