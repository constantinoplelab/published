function ProcessData_Figure6(datadir, codedir, savedir)
%ProcessData_Figure6 - Process raw data saved under datadir such that it can be plotted by PlotFigure6.
% INPUTS:
%   datadir - Local directory where 'Golden_Nature/RawData_Figure6.mat' from Zenodo was saved
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
load([datadir, 'RawData_Figure6'],...
    'TABLE_tetshEsr1', 'TABLE_lentivirus', 'lentirats', 'DoxBeh',...
    'DoxBeh_noscreen', 'ratTrial', 'SerumOsmolalityTable',...
    'EstradiolExp1Table', 'EstradiolExp2Table')

%% Process data
%Set general variables
ratlist_shEsr1 = {'G091', 'G093', 'G096', 'G097', 'G098',...
    'G133', 'G134', 'G135', 'G136', 'G139', 'G141'}; %G133 only has pre-virus cycling data

%--------------------------------------------------------------------------
% 6c. Quantification of RNAscope for tet-shEsr1
%--------------------------------------------------------------------------
TABLE_tet = TABLE_tetshEsr1;

%--------------------------------------------------------------------------
% 6d. Quantification of RNAscope for lenti-shEsr1
%--------------------------------------------------------------------------
TABLE_lenti = TABLE_lentivirus;

%--------------------------------------------------------------------------
% 6e. Example rat behavioral sensitivity to the blocks
%--------------------------------------------------------------------------
ITIbyBlock_shEsr1_screen = ITIbyBlock_dox_stages(ratlist_shEsr1,...
    lentirats, DoxBeh);

%--------------------------------------------------------------------------
% 6g. Model simulation of reduced gain on behavioral sensitivity to the
% reward blocks using G025's data 
%--------------------------------------------------------------------------
nexttile
Ls = cell(2,1);
for jj=1:2
    if jj==1 %pre
        params = [0.6 1]; %alpha, RPE gain factor
        gain_arg = 1;
    elseif jj==2 %shEsr1
        params = [0.6 0.25]; %alpha, RPE gain factor
        gain_arg=1;
    end
    %manipulating RPE gain
    Ls{jj} = GenerateLatencyData_VanillaAlpha_cg(params, ratTrial,... 
            0, gain_arg, 3);
end

%--------------------------------------------------------------------------
% 6h. Correlation of osmolality and estradiol in serum
%--------------------------------------------------------------------------
outlier_thresh = 2.5; %2.5*std outlier
[estradiol, RatConc] =...
    serum_osmolality_estradiol(outlier_thresh, SerumOsmolalityTable,...
    EstradiolExp1Table, EstradiolExp2Table); 

%--------------------------------------------------------------------------
% 6i. Population level shEsr1 affect on volume consumed
%--------------------------------------------------------------------------
ITIbyBlock_shEsr1_noscreen = ITIbyBlock_dox_stages(ratlist_shEsr1,...
    lentirats, DoxBeh_noscreen); %for volume consumed per session, we don't care about performance so don't screen


%% Save processed data
save([savedir 'ProcessData_Figure6'], 'ratlist_shEsr1',...
    'TABLE_tet', 'TABLE_lenti', 'ITIbyBlock_shEsr1_screen',...
    'ITIbyBlock_shEsr1_noscreen', 'ratTrial', 'Ls', 'estradiol',...
    'RatConc');

end