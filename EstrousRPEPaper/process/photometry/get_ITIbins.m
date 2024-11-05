function [ITIacrossrats, ITIbins_density] =...
    get_ITIbins(ratlist, Bstruct, numbins)    
%Use mixed blocks and violations only
%input
    %ratlist
    %numbins: number of bins to sort all ITIs into
 
ITIs_rats = cell(length(ratlist), 1);
numbins_vector = linspace(1,numbins,numbins);

for rat = 1:length(ratlist)    

    ratname = ratlist{rat};
    
    %Load photometry data     
    bdata = Bstruct.(ratname).CPIn;
    
    ITIs = bdata.ITI;
    Blocks = bdata.Block;        

    %remove ITIs > 98th percentile
    ITIs(ITIs>prctile(ITIs,98)) = NaN;    

    ITIs_rats{rat} = ITIs(Blocks==1 & bdata.PrevTrialType==2); %only use trials post violations
    
end

ITIacrossrats = cell2mat(ITIs_rats);
fractions = numbins_vector/numbins;

%% get ITI bins based on density
figure;
ITICDFplot = cdfplot(ITIacrossrats); hold on
Yvalues = ITICDFplot.YData;
Xvalues = ITICDFplot.XData;
ITIbins_density = NaN(length(fractions), 1);
for bin = 1:length(fractions)
    thisfract = fractions(bin);
    [~ , yfractidx] = min(abs(Yvalues-thisfract));
    Xatfract = Xvalues(yfractidx);
    xline(Xatfract, '--r');
    ITIbins_density(bin, 1) = Xatfract;
end
close all

end
