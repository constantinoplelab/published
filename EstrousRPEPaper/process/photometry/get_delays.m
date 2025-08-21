function [delaybins, delays_rats] = get_delays(ratlist, Bstruct, numbins)    
%run to find delay bins 

delays_rats = cell(length(ratlist), 1);
for rat = 1:length(ratlist)    
    ratname = ratlist{rat};
    bdata = Bstruct.(ratname).SideOff;
    thesestages = cellfun(@(x) logical(sum(strcmp(x,...
        {'Proestrus', 'Diestrus'}))), bdata.stage);
    delays = bdata.RewardDelay(thesestages);
    %remove catch trials
    noncatch_delays = delays;
    noncatch_delays(noncatch_delays==100) = NaN;
    delays_rats{rat} = noncatch_delays;
end

delays = [];
for rat = 1:length(ratlist)    
    delays = [delays; delays_rats{rat}];
end

numbins_vector = linspace(1,numbins,numbins);
fractions = numbins_vector/numbins;

%% get delay bins
figure;
delayCDFplot = cdfplot(delays); hold on
Yvalues = delayCDFplot.YData;
Xvalues = delayCDFplot.XData;
delaybins = NaN(length(fractions), 1);
for bin = 1:length(fractions)
    thisfract = fractions(bin);
    [~ , yfractidx] = min(abs(Yvalues-thisfract));
    Xatfract = Xvalues(yfractidx);
    xline(Xatfract, '--r');
    delaybins(bin, 1) = Xatfract;
end
close all
end
