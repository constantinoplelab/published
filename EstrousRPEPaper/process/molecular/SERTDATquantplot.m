function [Stages, AreaByStage, MeanFByStage, pvals] =...
    SERTDATquantplot(TABLE, protein)

if strcmp(protein, 'SERT')
    Stages = {'Proestrus', 'Estrus', 'Diestrus'};
    ratlist = {'RatA','RatB','RatC','RatD','RatE','RatF','RatG','RatH','RatI'};
    stages = {'Proestrus','Estrus','Diestrus','Proestrus','Estrus',...
        'Diestrus','Proestrus','Estrus','Diestrus'}; %identity for each rat in ratlist e.g., {'Proestrus', %'Estrus', 'Diestrus'}
    ch = 'ch00';
elseif  strcmp(protein, 'DAT')
    Stages = {'Proestrus', 'Diestrus'};
    ratlist = {'RatJ', 'RatK', 'RatL','RatM', 'RatN', 'RatO'};
    stages = {'Proestrus','Diestrus','Proestrus', 'Diestrus',...
        'Proestrus','Diestrus'}; %identity for each rat in ratlist e.g., {'Proestrus', %'Estrus', 'Diestrus'}
    ch = 'ch00';
elseif  strcmp(protein, 'TH')
    Stages = {'Proestrus', 'Estrus', 'Diestrus'};
    ratlist = {'RatA','RatB','RatC','RatD','RatE','RatF','RatG','RatH',...
        'RatI', 'RatJ', 'RatK', 'RatL','RatM', 'RatN', 'RatO'};
    stages = {'Proestrus','Estrus','Diestrus','Proestrus','Estrus',...
        'Diestrus','Proestrus','Estrus','Diestrus','Proestrus','Diestrus',...
        'Proestrus', 'Diestrus','Proestrus','Diestrus'}; %identity for each rat in ratlist e.g., {'Proestrus', %'Estrus', 'Diestrus'}
    ch = 'ch01';
end

%Remove RatG (outlier) 
TABLE{ismember(TABLE.RatID, 'RatG'), 7} = NaN;
TABLE{ismember(TABLE.RatID, 'RatG'), 8} = NaN;

numimages = 54; %the max number of images, use to pad with NaNs

AreaByStage = NaN(numimages, length(Stages));
MeanFByStage = NaN(numimages, length(Stages));
for s = 1:length(Stages)
    ratsthisstage = ratlist(ismember(stages, Stages{s}));
    
    %AREA
    areadata = TABLE.Area(ismember(TABLE.RatID,...
        ratsthisstage) & strcmp(TABLE.Fluoro, ch));    
    AreaByStage(1:length(areadata), s) = areadata; %rows are images, columns are stages

    %MEAN F
    Fdata = TABLE.MeanF(ismember(TABLE.RatID,...
        ratsthisstage));
    MeanFByStage(1:length(Fdata), s) = Fdata;

end
if ismember(protein, {'SERT','TH'})
    pvals.area.maineffect = kruskalwallis(AreaByStage,[], 'off');
    pvals.area.ProestrusEstrus = ranksum(AreaByStage(:,1),AreaByStage(:,2));
    pvals.area.ProestrusDiestrus = ranksum(AreaByStage(:,1),AreaByStage(:,3));
    pvals.area.EstrusDiestrus = ranksum(AreaByStage(:,2),AreaByStage(:,3));

    pvals.Fl.maineffect = kruskalwallis(MeanFByStage,[], 'off');
    pvals.Fl.ProestrusEstrus = ranksum(MeanFByStage(:,1),MeanFByStage(:,2));
    pvals.Fl.ProestrusDiestrus = ranksum(MeanFByStage(:,1),MeanFByStage(:,3));
    pvals.Fl.EstrusDiestrus = ranksum(MeanFByStage(:,2),MeanFByStage(:,3));

elseif ismember(protein, {'DAT'})

    pvals.area.maineffect = ranksum(AreaByStage(:,1),AreaByStage(:,2));
    pvals.Fl.maineffect = ranksum(MeanFByStage(:,1),MeanFByStage(:,2));

end

end