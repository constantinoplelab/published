function [RPEbins_density, RPEbins_equallyspaced] =...
    get_RPEbins(RPEs_rats, numbins)    
%Use mixed blocks and violations only
%input
    %numbins: number of bins to sort all RPEs into
 
numbins_vector = linspace(1,numbins,numbins);
RPEacrossrats = cell2mat(RPEs_rats);
fractions = numbins_vector/numbins;

%% get RPE bins based on density
f1 = figure;
ITICDFplot = cdfplot(RPEacrossrats); hold on
Yvalues = ITICDFplot.YData;
Xvalues = ITICDFplot.XData;
RPEbins_density = NaN(length(fractions), 1);
for bin = 1:length(fractions)
    thisfract = fractions(bin);
    [~ , yfractidx] = min(abs(Yvalues-thisfract));
    Xatfract = Xvalues(yfractidx);
    xline(Xatfract, '--r');
    RPEbins_density(bin, 1) = Xatfract;
end
close(f1)

RPEbins_equallyspaced = linspace(-3, 3, numbins+1);

figure;
cdfplot(RPEacrossrats); hold on
for bin = 1:length(RPEbins_equallyspaced)
    xline(RPEbins_equallyspaced(bin), '--r'); hold on
end

close all

end
