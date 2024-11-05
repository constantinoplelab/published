function [avgITIbyblock_opto, errITIbyblock_opto,...
    avgITIbyblock_nonopto, errITIbyblock_nonopto] =...
    CompareOptoITIbyBlock(optoafterdate, Aligned, S) %controlafterdate, varargin
%input: ratlist (e.g. {'G060'})
        %afterdate is the date after opto surgery for control sessions 
            %G071: optoafterdate: 20221207, control: '2022-11-06', remove
            %one session after control?
            %G072: optoafterdate: 20221207, control: '2022-11-06'
            %G073: optoafterdate=20230214 for 10 mW,  controlafterdate='2022-10-18' (check for CPIn opto)
            %G081: optoafterdate=20230305, controlafterdate='2023-02-22'
            %G087: optoafterdate=20230205, controlafterdate='2023-01-23',
            %after 1 week after recovering from surgery
            %G089: optoafterdate=20230212 (20230314 for 7 mW), controlafterdate='2023-01-30'
            %after 1 week after recovering from surgery
            %G095: optoafterdate=20230601, controlafterdate='2023-05-23'
            %after 1 week after recovering from surgery
            %G102: optoafterdate=20230523, controlafterdate='2023-05-17' 
            %after 1 week after recovering from surgery
            %G103: optoafterdate=20230523, controlafterdate='2023-05-17' 
            %after 1 week after recovering from surgery
            %G102 and G103: optoafterdate=20230523, controlafterdate='2023-05-17',
            %G095: optoafterdate=20230601, controlafterdate='2023-05-23',
            %after 1 week after recovering from surgery,first day getting
            %plugged in
            %G092: optoafterdate=20230705, controlafterdate='2023-06-22',
            %after 1 week after recovering from surgery,first day getting
            %plugged in

prctile_cutoff = 98;

%create functions
block_avg = @(x, blocks) arrayfun(@(y) mean(x(blocks==y), 'omitnan'), 1:3);
block_err = @(x, blocks) arrayfun(@(y)...
    std(x(blocks==y), 'omitnan')./sqrt(sum(~isnan(x(blocks==y)))), 1:3);

%% OPTO data
%use dates based on when opto started
datesnums = cellfun(@str2num,Aligned.SessDate);
optoafterdatenum = str2num(optoafterdate);
Aligned = Aligned(datesnums > optoafterdatenum, :);

%remove opto sessions with SideOff stimulation (safety precaution)
SideOff_dates = unique(Aligned.SessDate(strcmp(Aligned.OptoEvent, 'SideOff')));
Aligned = Aligned(~ismember(Aligned.SessDate, SideOff_dates), :);
ITIs = Aligned.ITI;

%Find number of trials per session
uniquedates = unique(Aligned.SessDate);
ntrials = NaN(length(uniquedates),1);
for sess = 1:length(uniquedates)
    ntrials(sess) = sum(strcmp(Aligned.SessDate, uniquedates(sess))); %find number of each session
end

%Process ITIs
%%number for sessions from this stage
ITIs(cumsum([0; ntrials(1:end-1)])+1) = NaN;
ITIs(ITIs>prctile(ITIs, prctile_cutoff)) = NaN;
%convert number of stimulations to logical
%split ITIs into sessions and remove effect of satiety over session
counter = 0;
lat_cell = cell(length(ntrials), 1);
block_cell = cell(length(ntrials), 1);
for j = 1:length(ntrials)
    lat_cell{j} = ITIs(counter+1:counter+ntrials(j));
    block_cell{j} = Aligned.Block(counter+1:counter+ntrials(j));
    counter = counter + ntrials(j);
end
ITIs_opto_processed = cell2mat(lat_cell);
Blocks_opto = cell2mat(block_cell);

%% CONTROL data
is_struct = cell2mat(cellfun(@isstruct,S.pd,'UniformOutput',false));
goods = cellfun(@(s) sum(s.TrainingStage==9)==length(s.TrainingStage), ...
    S.pd);
S.pd = S.pd(is_struct & goods);
S.peh = S.peh(is_struct & goods);
numnonoptosessions = length(S.pd);
%calculate ITI vector across sessions and NaN first trial
ITIs_over_sessions_nonopto = getITI(S); %includes post-violations
for j = 1:numnonoptosessions
    ITIs_thissess = ITIs_over_sessions_nonopto{j};
    ITIs_NaNfirsttrial = [NaN; ITIs_thissess(2:end)];
    ITIs_over_sessions_nonopto{j} = ITIs_NaNfirsttrial;
end
ITIs = cell2mat(ITIs_over_sessions_nonopto');
%Find ITI threshold to cutoff
lat_thresh = prctile(ITIs,prctile_cutoff);
lat_cell = cell(numnonoptosessions,1);
blocks_over_sessions_nonopto = cell(numnonoptosessions, 1);
%remove outliers and effect of satiety over session
for j = 1:numnonoptosessions
    ITIs_processed = ITIs_over_sessions_nonopto{j};
    ITIs_processed(ITIs_processed>lat_thresh) = NaN;
    lat_cell{j} = ITIs_processed;
    blocks_over_sessions_nonopto{j} = S.pd{j}.Block;
end
ITIs_nonopto_processed = cell2mat(lat_cell);
Blocks_nonopto = cell2mat(blocks_over_sessions_nonopto);

%calculate averages
avgITIbyblock_opto = block_avg(ITIs_opto_processed, Blocks_opto);
errITIbyblock_opto = block_err(ITIs_opto_processed, Blocks_opto);

avgITIbyblock_nonopto = block_avg(ITIs_nonopto_processed, Blocks_nonopto);
errITIbyblock_nonopto = block_err(ITIs_nonopto_processed, Blocks_nonopto);

end