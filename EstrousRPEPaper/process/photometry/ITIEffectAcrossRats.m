function [da_rats, AUC_overrats, P_byrat, Rcorr_byrat, bins, T,...
    ITI_legend] = ITIEffectAcrossRats(ratlist, ITIbins_density, binnum,...
    event, window, Pstruct, Bstruct)
%Compares reward response to offer cue BEFORE different initiation times
%input:
    %binnum = number of bins to split initiation times into

da_rats = cell(length(ratlist), 1);
AUC_overrats = NaN(length(ratlist), binnum);
P_byrat = NaN(length(ratlist),1);
Rcorr_byrat = NaN(length(ratlist),1);
for rat = 1:length(ratlist)    
    
    ratname = ratlist{rat};
    disp(['Initiation time photometry: ' ratname])

    [DA, AUC_overrats(rat,:), P_byrat(rat,1), Rcorr_byrat(rat,1), T,...
        ITI_legend, bins] = ITI_effect(Bstruct.(ratname),...
        Pstruct.(ratname), ITIbins_density, event, binnum, window);
            
    da_rats{rat} = DA;

end

end
