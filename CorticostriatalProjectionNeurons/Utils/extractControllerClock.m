function [controller_start, controller_fs] = extractControllerClock(fileID)
% Collects the start and sample rate of the recording stream from the saved
% test file sync_messages. Updated 08/2024 SSS

% Input: 
%   fileID: file to open

% Output:
%   controller_start: time of data stream start in sample number
%   controller_fs: intended sample rate of the data stream

%% Open sync_messages.txt
cstring = fscanf(fileID,'%s');
% if contains(cstring,'Softwaretime')
if contains(cstring, 'Processor')
    %% OpenEphys version 0.5
    nstring = split(cstring,'Hz'); %end chars in sync_messages
    a = 'starttime:';
    b = '@';
    start = 1; %controller start is first
    rate = 2; %sample rate is second
else
    nstring = split(cstring,'StartTimefor'); %start char in sync_messages
    a = '@';
    b = 'Hz:';
    start = 2; %controller start is second
    rate = 1; %sample rate is first
end

for j = 1:numel(nstring)

    %extract controller clock
    if strfind(nstring{j}, 'NI-DAQ') >0
        cj = nstring{j};
        sstart = strfind(cj, a);
        data = cj(sstart+numel(a):end);
        data_split = split(data, b);
        controller_start.daq = str2double(data_split{start});
        controller_fs.daq = str2double(data_split{rate});

        %extract probe clocks -- as of 8/12/2024 these can't be changed
        %anyways and should always be 30k for the AP stream, 2.5k for the
        %LFP stream
    elseif strfind(nstring{j}, 'Neuropix-PXI')>0
        if strfind(nstring{j},'2500')>0 %by definition this will have to be 2500 so not really necessary
            cj = nstring{j};
            sstart = strfind(cj, a);
            data = cj(sstart+numel(a):end);
            data_split = split(data, b);
            controller_start.lfp = str2double(data_split{start});
            controller_fs.lfp = str2double(data_split{rate});
        else
            cj = nstring{j};
            sstart = strfind(cj, a);
            data = cj(sstart+numel(a):end);
            data_split = split(data, b);
            controller_start.ap = str2double(data_split{start});
            controller_fs.ap = str2double(data_split{rate});
        end
    end

end