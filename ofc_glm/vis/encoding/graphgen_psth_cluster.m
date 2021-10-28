function [f,f1] = graphgen_psth_cluster(psth,clustlabs,clustfilename,useZscore,varargin)
%simple script for calculating PSTH or new clusteirng setups

%useZscore = true;

if numel(varargin) > 0
    ops = varargin{1};
    useSEM = ops.SEM;
    if isfield(ops,'fignum')
        fignum = ops.fignum;
    else
        fignum = 1;
    end
else
    ops = struct('vishandle','true');
    useSEM = true;
    fignum = 1;
end

%define a separate clustering order, to average over one labeling,
%but sort in another
if isfield(ops,'clustlabs_sort')
    clustlabs_sort = ops.clustlabs_sort;
else
    clustlabs_sort = clustlabs;
end

if useZscore
    %load psth_dat first from an anlauysis fpolder
    %psth = DatPSTH_dat{1}(:,:,1);
    meandat = nanmean(psth,2);
    stddat = nanstd(psth,[],2);
    psth = (psth-meandat)./stddat;
    ytext = 'rate (z-scored)';
else
    %load psth_dat first from an anlauysis fpolder
    %psth = DatPSTH_dat{1}(:,:,1);
    ytext = 'rate (Hz)';
end

f = figure();
clf
nclust = numel(unique(clustlabs));
%dirty hack for right now. 
if size(psth,2) ==  121
    tmesh = -2:0.05:4;
elseif size(psth,2) ==  161   
    tmesh = -4:0.05:4;
else
    disp('this sizing of psth time domain not supported')
end


if numel(varargin) > 1
    colormat = varargin{2};
else
    if nclust < 10
        colormat = linspecer(nclust,'qualitative');
    else
        colormat = linspecer(nclust);
    end
end



for j = 1:nclust

    if useSEM
        
        if sum(clustlabs==j-1) >1
            sem = nanstd(psth(clustlabs==j-1,:))/sqrt(sum(clustlabs==j-1));
            
            shadedErrorBar(tmesh,nanmean(psth(clustlabs==j-1,:)),sem,...
            'lineprops',{'color',colormat(j,:),'linewidth',2})
        else
             plot(tmesh,psth(clustlabs==j-1,:),'color',colormat(j,:),'linewidth',2)
        end
    
    else
        if sum(clustlabs==j-1) >1
            plot(tmesh,nanmean(psth(clustlabs==j-1,:)),'color',colormat(j,:),'linewidth',2)
        else
            plot(tmesh,psth(clustlabs==j-1,:),'color',colormat(j,:),'linewidth',2)
        end
            
    end
        
    hold on
end

    
xlabel('time to start')
legend(num2str((1:nclust)'),'location','best')
ylabel(ytext);


namebase = split(clustfilename,'/');
namebase = namebase{end};
namebase = split(namebase,'.mat');
namebase = namebase{1};
if useZscore
    title(strcat('cluster PSTH: ',namebase,' zscore'),'interpreter','none')
else
    title(strcat('cluster PSTH: ',namebase),'interpreter','none')
end
set(gca,'fontsize',15)
vline(0,'k')


%plot the heat mat
f1 = figure()
clf
nn = size(psth,1);
[I,groupbreaks] = time2peakByCluster(psth,struct('grouping',clustlabs_sort));
imagesc(tmesh,1:nn,psth(I,:))
hold on

for k = 1:numel(groupbreaks)
    hline(groupbreaks(k),'r')
end

colorbar()
xlabel('time to start (s)')
ylabel('neurons')
title('PSTH (z-score)') 
set(gca,'fontsize',15)
vline(0,'k');




    