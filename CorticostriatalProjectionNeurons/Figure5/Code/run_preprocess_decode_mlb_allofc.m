%% the preprocessing of data for all-ofc MLB decoding. 
% so that actual fitting can run on torch 

%
%savedir = '/Users/dhocker/projects/dynamics/data/maggie/';
%savename = strcat(savedir,'preprocess_second_mlb_alldat_allOFC_quantile_psthsubtract.mat');

function [] = run_preprocess_decode_mlb_allofc(savename, epoch, decodertype)

    %location of analysis codebase
    codepath = '~/projects/constantinoplelab/Analysis/';
    
    % add code paths
    addpath(genpath(codepath))
    addpath(genpath(codepath + "david"))
    addpath(genpath(codepath + "maggie"))
    %% load table for finding sessions
    Etable = load('/Users/dhocker/projects/dynamics/data/maggie/EphysTable.mat');
    ETable = Etable.ETable;
    
    test = ETable{:,"recording_site"};
    mask1 = cellfun(@(x) strcmp(x,'OFC'),test,'UniformOutput',true);
    
    test = ETable{:,"fiber_site"};
    mask2 = cellfun(@(x) strcmp(x,'DLS'),test,'UniformOutput',true);
    %this does not matter
    
    test = ETable{:,"protocol"};
    mask3 = cellfun(@(x) strcmp(x,'RWTautowait2'),test,'UniformOutput',true);
    %mask3 = cellfun(@(x) strcmp(x,'RWTautowait2OptoTest'),test,'UniformOutput',true);
    % does not currwently matter
    
    test = ETable{:,"stimulation"};
    %mask4 = test==1; % stimulation
    mask4 = test==0; % no stimulatio

    goodsessions = find(mask1 & mask2 & mask3 & mask4);
    %goodsessions = find(mask1 & mask4);
    nsess = numel(goodsessions);
    
    disp('number of sessions')
    disp(nsess)

    %% get filenames for processed data of  good sessions
    
    loadmac = true;
    
    fullnames = cell(nsess,1);
    
    for j = 1:nsess
    
        if loadmac
            usedir = ETable{goodsessions(j),'savepath'}{1};
            usedir = replace(usedir,'\\constantinoplelab.cns.nyu.edu\server2\','/Volumes/server2/');
            usedir = replace(usedir,'\\constantinoplelab.cns.nyu.edu\server3\','/Volumes/server3/');
            usedir = replace(usedir,'\','/');
    
        else
            usedir = ETable{goodsessions(j),'savepath'}{1};
        end
    
        
        fullnames{j}  = strcat(usedir,ETable{goodsessions(j),'matfile'}{1});
        disp(ETable{goodsessions(j),'matfile'}{1})
    end
    
    % if .mat file was not present, remove from cell and get new nsess
    goodsess = cellfun(@(x) numel(strfind(x,".mat")) > 0,fullnames);
    fullnames = fullnames(goodsess);
    nsess = numel(fullnames);

%%  get session data, count number neurons
    disp('retrieving data for each session')
    %Sall = cell(nsess,1);
    %for j = 42:nsess
    for j = 1:nsess
    %for j = 1:1
        disp(j)
        disp(fullnames{j})
        sj = struct();
        [output] = decode_parsedata_rawdat(fullnames{j},epoch, decodertype);
        sj.output = output;
        sj.id = j;
    
        %check if there are enough samples per catetory: 
        ntrials_percond = 20; %number of trials per condition
        mlb_m = sum(sj.output.mostlikely_block == 1 & sj.output.rewarded_trials'==1);
        mlb_h = sum(sj.output.mostlikely_block == 2 & sj.output.rewarded_trials'==1);
        mlb_l = sum(sj.output.mostlikely_block == 3 & sj.output.rewarded_trials'==1);
        mlb2_m = sum(sj.output.secondlikely_block == 1 & sj.output.rewarded_trials'==1);
        mlb2_h = sum(sj.output.secondlikely_block == 2 & sj.output.rewarded_trials'==1);
        mlb2_l = sum(sj.output.secondlikely_block == 3 & sj.output.rewarded_trials'==1);


            if mlb_m >= ntrials_percond && mlb_h >= ntrials_percond && mlb_l >= ntrials_percond
                % check mlb2 as well
                if mlb2_m >= ntrials_percond && mlb2_h >= ntrials_percond && mlb2_l >= ntrials_percond
                    Sall{j} = sj;
                else
                    Sall{j} = [];
                end
            else
                Sall{j} = [];
            end
    end
    
    %remove bad sessions
    badsess = cellfun(@(x) size(x,1),Sall,'UniformOutput',true) == 0;
    Sall(badsess) = [];
    %nsess = numel(Sall);

    %% save
    disp('saving')
    save(savename,'*');

end


