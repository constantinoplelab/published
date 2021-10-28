function [condmask,condnames] = condmaskLibrary(type,file_idx,A,ops)
%create a conditional mask library for different filtering operations
%
%Inputs:
%   type: string to decide what to condition on
%   file_idx: neuron index in concatenated A data structure, or index in
    %           fnames
%   A: concatenated alldata structure that holds all of the ofc data
%   ops: extra options structure in case conditional requires a numerical
%           range, such as reward volume
%
%Output:
%   condmask: boolean of which trials in A{-}.trial_idx have that condition


%find correct neuron
for j = 1:numel(A)
    if A{j}.fileind==file_idx
        Dat = A{j};
        break
    end
end


%some preliminaries that are common for a few of the cases
water = (Dat.hits).*Dat.chosenval(1:length(Dat.hits)); %water received on each trial
water(isnan(water))=0; %omit nans from opt outs. change to no water
water_prev = [0;water(1:end-1)];
rewRate = cumsum(water,'omitnan')./(1:numel(water))';
rewRate(isnan(rewRate))=0; %some nans at the beginning
prevRewRate = [0;rewRate(1:end-1)]; %average reward up until previous trial

%prospect theory form
chosenprob = Dat.chosenprob;
chosenprob(isnan(chosenprob)) = 0;
RPE = water-water.*chosenprob; %reward prediction error on current trial
prevRPE = [0;RPE(1:end-1)]; %RPE of previous trial

%session progress
session_progress = (1:numel(Dat.handles.start))/(numel(Dat.handles.start));


%find which condition you want
switch type
    case 'all'
        condmask = true(size(Dat.trial_idx));
        condnames = type;
    case 'prevWin'
        condmask = Dat.prev_hits==1;
        condnames = type;
    case 'prevLoss'
        condmask = Dat.prev_hits==0;
        condnames = type;
    case 'prevOptOut'
        condmask = isnan(Dat.prev_hits);
        condnames = type;
    case 'left'
        condmask = Dat.went_right==0;
        condnames = type;
    case 'right'
        condmask = Dat.went_right==1;
        condnames = type;
    case 'win'
        condmask = Dat.hits==1;
        condnames = type;
    case 'loss'
        condmask = Dat.hits==0;
        condnames = type;
    case 'risky'
        isRightSafe=  Dat.right_prob==1;
        isLeftSafe=   Dat.left_prob==1;    
        choseSafe = (isRightSafe & Dat.went_right==1) | ...
                (isLeftSafe & Dat.went_right==0);
        %choseRisky = (isRightSafe & Dat.went_right==0) | ...
                %(isLeftSafe & Dat.went_right==1);
            
        choseRisky = ~choseSafe & ~isnan(Dat.hits);
            
        condmask = choseRisky ==1;
        condnames = type;
    case 'safe'
        isRightSafe=  Dat.right_prob==1;
        isLeftSafe=   Dat.left_prob==1;    
        choseSafe = (isRightSafe & Dat.went_right==1) | ...
                (isLeftSafe & Dat.went_right==0);
            
        choseSafe = choseSafe & ~isnan(Dat.hits);
            
         condmask = choseSafe == 1;
         condnames = type;
    case 'vol'
        %discretize into M equal parts and mask.
        %how to save masks though? as cell array 
        condbins = [0,6,12,24,48];
        condnames = {'0 ul','6 ul','12 ul','24 ul','48 ul'};
        ncond = numel(condbins);
        condmask = cell(ncond,1);
        for j = 1:ncond
            condmask{j} = water == condbins(j);
        end
        
    case 'prob'
        %TODO: CODE up chosen probability
       
    case 'prevVol'
        condbins = [0,6,12,24,48];
        condnames = {'0 ul','6 ul','12 ul','24 ul','48 ul'};
        ncond = numel(condbins);
        condmask = cell(ncond,1);
        for j = 1:ncond
            condmask{j} = water_prev == condbins(j);
        end
        
    case 'rewardRate'
        ncond = 5;
        condbins = linspace(min(rewRate),max(rewRate)+1,ncond+1); %+1 fixes edge case of inequality test below
        condnames = cellfun(@(x) num2str(x),num2cell(condbins),'UniformOutput',false);
        condmask = cell(ncond,1);
        for j = 1:ncond
            condmask{j} = (rewRate >= condbins(j)) & (rewRate < condbins(j+1));
        end
        
    case 'prevRewardRate'
        ncond = 5;
        condbins = linspace(min(prevRewRate),max(prevRewRate),ncond+1); %+1 fixes edge case of inequality test below
        condnames = cellfun(@(x) num2str(x),num2cell(condbins),'UniformOutput',false);
        condmask = cell(ncond,1);
        for j = 1:ncond
            condmask{j} = (prevRewRate >= condbins(j)) & (prevRewRate < condbins(j+1));
        end
        
     case 'RPE'
        ncond = 5;
        condbins = linspace(min(RPE),max(RPE)+1,ncond+1); %+1 fixes edge case of inequality test below
        condnames = cellfun(@(x) num2str(x),num2cell(condbins),'UniformOutput',false);
        condmask = cell(ncond,1);
        for j = 1:ncond
            condmask{j} = (RPE >= condbins(j)) & (RPE < condbins(j+1));
        end   
        
    case 'prevRPE'
        ncond = 5;
        condbins = linspace(min(prevRPE),max(prevRPE)+1,ncond+1); %+1 fixes edge case of inequality test below
        condnames = cellfun(@(x) num2str(x),num2cell(condbins),'UniformOutput',false);
        condmask = cell(ncond,1);
        for j = 1:ncond
            condmask{j} = (prevRPE >= condbins(j)) & (prevRPE < condbins(j+1));
        end 
        
    case 'sessProg'
        ncond = 5;
        condbins = linspace(min(session_progress),1.01,ncond+1); %+1 fixes edge case of inequality test below
        condnames = cellfun(@(x) num2str(x),num2cell(condbins),'UniformOutput',false);
        condmask = cell(ncond,1);
        for j = 1:ncond
            condmask{j} = (session_progress >= condbins(j)) & (session_progress < condbins(j+1));
        end  
        
    case 'vol_no0'
        %discretize into M equal parts and mask.
        %how to save masks though? as cell array 
        condbins = [6,12,24,48];
        condnames = {'6 ul','12 ul','24 ul','48 ul'};
        ncond = numel(condbins);
        condmask = cell(ncond,1);
        for j = 1:ncond
            condmask{j} = water == condbins(j);
        end
        
    case 'prevVol_no0'
        condbins = [6,12,24,48];
        condnames = {'6 ul','12 ul','24 ul','48 ul'};
        ncond = numel(condbins);
        condmask = cell(ncond,1);
        for j = 1:ncond
            condmask{j} = water_prev == condbins(j);
        end
        
    otherwise %assume no masking
        disp(strcat(['not found: ', type]))
        condmask = true(size(Dat.trial_idx));
        
end
  
if isfield(ops,'removeOptOut')
    if ops.removeOptOut
        if iscell(condmask)
            for j = 1:numel(condmask)
                condmask{j} = condmask{j}(~isnan(Dat.hits));
            end
        else
            condmask = condmask(~isnan(Dat.hits));
        end
    end
end
       
        
        