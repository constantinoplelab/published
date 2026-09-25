function ops = make_ops(ETable)
% make ops struct to store info about raw data
% INPUTS: shared ephys table with info for each session, only containing
% the session of interest
% OUTPUTS: ops struct with information about # of channels, sampling rate,
% index of channels used for syncing, opto LED times etc. 
% dh: Dec. 2019. SS updated 2024

APpath = ETable.fullpath{1};
LFPpath = ETable.LFPpath{1};
ops = struct();

ops.ratname = ETable.ratname;
ops.sessiondate = ETable.sessiondate;
ops.lfp = LFPpath;
ops.ksdir = APpath;
name = split(APpath, '\');
ops.datadir = strjoin(name(1:end-2),filesep);
fpath = strjoin(name(1:end-2),filesep);

%ops.datadir = fpath; %main head directory for neuropixels data
%ops.ksdir = fullfile(ops.datadir,'Record Node 103/experiment1/recording1/continuous/Neuropix-PXI-100.0'); %where your kilosort data lives
%fpath = fullfile(fpath, 'Record Node 103/experiment1/recording1');

ops.nCh = 384; %number of recorded channels from probe
%ops.fs = 30000; %sampling frequency of probe, in Hz (MD removed to select
%from sync messages instead of being hard coded. 
ops.binfile = 'continuous.dat'; %where waveforms live in ops.ksdir

%sync information
%1-indexed digital channel input on controller to poke ports 
%(note: most controllers used zero-indexing, and you need to add 1)
ops.Ch_timesync = 1;
ops.Ch_cpoke = 2; 
ops.Ch_lpoke = 6;
ops.Ch_rpoke = 3;
ops.Ch_sLED = 8;    %right for optogenetic stimulation
ops.Ch_sLEDL = 7;   % left for optogenetic stimulation NOT optotagging
ops.Ch_trialnumber = 1; %trial number for analog trial number. 1-indexed
ops.controllerName = 'NI-DAQmx'; %name of contorller sendin sync signals in sync_messages.txt
ops.probeName = 'Neuropix-PXI'; %similar name of probe
ops.probeProc_clock = 0; %subprocess within that probe. 0 should be AP, 1 should be LFP.

if strcmp(ETable.ratname, 'G160') & ETable.sessiondate==datetime('2026-04-15') | ETable.sessiondate==datetime('2026-04-16') %bncs for center and right were accidentally swapped
    ops.Ch_cpoke = 3; 
    ops.Ch_rpoke = 2;
end

%find event clock start time (saved in sync_messages)
fileID = fopen(fullfile(fpath,'sync_messages.txt'),'r');    % changed with my settings md

[~, controller_fs] = extractControllerClock(fileID);
ops.controller_fs = controller_fs.daq; %NI-DAQ sample rate
if isfield(controller_fs, 'ap')
    ops.probe_fs = controller_fs.ap; %Probe sample rate, AP stream
else
    ops.probe_fs = 30000; %only NI-DAQ info was saved in sync messages for earlier versions, but sample rate is always 30k
end

if isfield(controller_fs, 'lfp')
    ops.lfp_fs = controller_fs.lfp; %Probe sample rate, LFP stream
else
    warning('No separate LFP channel found; sampling rate set to probe sampling rate.')
    ops.lfp_fs = ops.probe_fs; 
    % Revised AP 4/22/26 - We think that it's wrong to have a hard-coded
    %   sampling rate of 2500. Hopefully won't cause issues
    %ops.lfp_fs = 2500; %only NI-DAQ info was saved in sync messages for earlier versions, but sample rate is always 2.5k
end

%waveform analysis parameters. not needed right now
ops.numWFperSess = 2; %number of waveforms persession to show. to assess drift
ops.numWFav = 2000; % number of spikes to average for each waveform 
ops.numchannels = 1; %number of channels to show the average waveforms on. 0 shows all
ops.wfWin = [-40,41]; % # samples to keep before and after spike to generate waveform. 
%default is [-40,41], whih for 30k sampling is 2.7 ms   


try 
    npx_dir = 'C:\Users\mld9131\Documents\GitHub\npx_analysis'; %where your npx_anlaysis is %%CHANGED THIS
    %add files to repo                    
    addpath(genpath(fullfile(npx_dir,'repos','spikes')));
    addpath(genpath(fullfile(npx_dir,'repos','npy-matlab')));
    addpath(genpath(fullfile(npx_dir,'repos','analysis-tools-openephys')));
    addpath(genpath(fullfile(npx_dir,'utils')));
catch
    disp("where is your npx_analysis at?")
end                      



end
