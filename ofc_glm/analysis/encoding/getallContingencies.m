function [C,volLR,probLR,condlist] = getallContingencies(file_idx,A,ops)
  
%gets contingencies for all files from A struct
condlist = {'prevWin'...,
                'prevLoss',...
                'prevOptOut',...
                'left',...
                'right',...
                'win',...
                'loss',...
                'risky',...
                'safe'};
%'rewHist',...
%'rewVol',...
%'RPE'};

   
%do a first one to get sizing
ctest = condmaskLibrary('win',file_idx,A,ops);
ns = numel(ctest);                
ncond = numel(condlist);

C = zeros(ncond,ns);

for j = 1:ncond    
    C(j,:) = condmaskLibrary(condlist{j},file_idx,A,ops);
end

%get water and prob on all trials
for j = 1:numel(A)
    if A{j}.fileind==file_idx
        Dat = A{j};
        break
    end
end

volLR = [Dat.this_left_volume,Dat.this_right_volume];
probLR = [Dat.left_prob, Dat.right_prob];

if isfield(ops,'removeOptOut')
    if ops.removeOptOut
        volLR = volLR(~isnan(Dat.hits),:);
        probLR = probLR(~isnan(Dat.hits),:);
    end
end



