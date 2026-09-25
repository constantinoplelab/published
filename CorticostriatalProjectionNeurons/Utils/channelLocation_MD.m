function SU = channelLocation_MD(S,SU)
%% Updates the SU to include the AP, ML, & DV of the electrode. 
% Note that this assumes that the electrode goes straight up without changing
% the AP/ML. Because there is definitely some angle associated with implants,
% I widened the range of acceptable values for the LO from what shannon had, 
% and chose reasonable(?) ranges for the DS(-0.2) and VS(1.2).
% inputs: 
%     S: S Struct for a single session
%     SU: SU cell array of every cell to analyse
%     tablepath: fullpath to your Npxl table
% outputs: 
%     SU.AP: AP relative to bregma of electrode given the table input of the shank
%     SU.ML: ML relative to bregma of the electrode given the table input of the shank
%     SU.DV: DV relative to bregma along the shank
%     SU.location: word description of the brain region the electrode is located
% At some point this should get updated to include adjustments for
% neuropixels 2.0, but it will need the implant orientation

load('\\constantinoplelab.cns.nyu.edu\server2\PhysiologyData\EphysTable.mat');
T = ETable;
try
isRat = find(strcmp(T.ratname,S.RatName) & T.sessiondate==S.SessionDate & strcmp(T.recording_site,S.recording_site) & T.session_num==S.session_num);
catch
    isRat = find(strcmp(T.ratname,S.RatName) & T.sessiondate==S.SessionDate & T.session_num==1);
end

%shank depth in bregma coordinates
DV_shank = T.DV_probe(isRat(1));
if isnan(DV_shank)
    disp('Session does not have a DV listed')
    return;
end

%cutoff for each region
if strcmp(T.recording_site{isRat},'OFC')
    cutoff1 = 6.1; %I've widened this range from Shannon
    cutoff2 = 3.7; %Again widened the range
    for i = 1:length(SU)
        DV_electrode = DV_shank-(SU{i}.channel_depth/1000);

        SU{i}.AP = T.AP_probe(isRat(1));
        SU{i}.ML = T.ML_probe(isRat(1));
        SU{i}.DV = DV_electrode;
try 
    DV_electrode<=cutoff1 && DV_electrode>=cutoff2;
catch
    keyboard
end

        if DV_electrode>cutoff1 %ventral to OFC
            SU{i}.location = 'Pir';
        elseif DV_electrode<=cutoff1 && DV_electrode>=cutoff2 %within OFC
            SU{i}.location = 'OFC';
        elseif DV_electrode<cutoff2 %dorsal to OFC
            SU{i}.location = 'white';
        end
    end

elseif strcmp(T.recording_site{isRat},'DLS')
    cutoff1 = 5.8; %
    cutoff2 = 3.4; %
    for i = 1:length(SU)
        DV_electrode = DV_shank-(SU{i}.channel_depth/1000);

        SU{i}.AP = T.AP_probe(isRat(1));
        SU{i}.ML = T.ML_probe(isRat(1));
        SU{i}.DV = DV_electrode;

        if DV_electrode>cutoff1 %ventral to DLS
            SU{i}.location = 'VS';
        elseif DV_electrode<=cutoff1 && DV_electrode>=cutoff2 %within DS
            SU{i}.location = 'DS';
        elseif DV_electrode<cutoff2 %dorsal to DLS
            SU{i}.location = 'Cortex';
        end
    end
elseif strcmp(T.recording_site{isRat},'VS')
    cutoff1 = 8; %
    cutoff2 = 5.8; %
    for i = 1:length(SU)
        DV_electrode = DV_shank-(SU{i}.channel_depth/1000);

        SU{i}.AP = T.AP_probe(isRat(1));
        SU{i}.ML = T.ML_probe(isRat(1));
        SU{i}.DV = DV_electrode;

        if DV_electrode>cutoff1 %ventral to VS
            SU{i}.location = 'Entorhinal';
        elseif DV_electrode<=cutoff1 && DV_electrode>=cutoff2 %VS
            SU{i}.location = 'VS';
        elseif DV_electrode<cutoff2 %dorsal to VS
            SU{i}.location = 'DS';
        end
    end
end
