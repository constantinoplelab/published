function  [EMdata, FracSIGIC, FracSIGMem, groups, groupID] =...
    DATlocation(EMdata)

%Get important variables
%make all ambiguous into membranous
EMdata.SIGCategory(strcmp(EMdata.SIGCategory, 'ambiguous')) = {'membranous'}; %DAT target (N-terminus) is intracellular

ratlist = unique(EMdata.RatID);
EMdata.GroupID = cell(size(EMdata, 1), 1);
for row = 1:size(EMdata, 1)
    if ismember(EMdata.RatID(row), {'EM4', 'EM6', 'EM7', 'EM9', 'EM10', 'EM11'})
        EMdata.GroupID(row) = {'estrus'};
    elseif ismember(EMdata.RatID(row), {'EM1', 'EM3', 'EM8', 'EM12', 'EM13', 'EM14'})
        EMdata.GroupID(row) = {'diestrus'};
    end
end

%Separate by category
FracSIGIC = NaN(length(ratlist), 1);
FracSIGMem = NaN(length(ratlist), 1);

%get stats for each rat
for rat = 1:length(ratlist)
    ratT = EMdata(strcmp(EMdata.RatID, ratlist(rat)), :);
    % SIG
    ratT.SIGMinimumDistancetoMembrane(strcmp(ratT.SIGCategory, 'ambiguous')) = 0; %DAT target (N-terminus) is intracellular    
    ratT.SIGCategory(strcmp(ratT.SIGCategory, 'ambiguous')) = {'membranous'}; %DAT target (N-terminus) is intracellular
    FracSIGIC(rat, 1) = sum(strcmp(ratT.SIGCategory, 'intracellular'))/size(ratT, 1);
    FracSIGMem(rat, 1) = sum(strcmp(ratT.SIGCategory, 'membranous'))/size(ratT, 1);
end

%%Plot by group
[~, idx] = unique(EMdata.RatID);
groupID = EMdata.GroupID(idx);
groups = flip(unique(groupID));

end