function [trials,Ftrials] = volbin_all_cells(all_SU,goods,all_S,all_index,win,varargin)
%plot heatmats aligned to reward volume for all cells 5/10 and 40/80 binned

% input:
% all_SU - the data of all cells
% all_behEvents - data of all cells
% all_index - indexing of all cells
% win - [start,end] of window you want plotted
% varargin:
    % alignto to sort the data too, otherwise defaults to time to peak ff on COff
    % plot or not

% output:
% figure of heatmat with z-scored ff in alinged to reward volume
% trials: low, mix, high
% ftrials: trials following the current volume
datpath = 'C:\Users\mld9131\Documents\GitHub\Papers\MLD\Figure2\Data';
load(fullfile(datpath,'TriColor.mat'))

%% z-score each cell, each trial
nanzscore = @(x)(x - mean(x,2, 'omitnan'))./std(x,0,2, 'omitnan');
alignto = {'COFF', 'SON','SOFF','Rew','Opt','CON'};
caxis = [-2,2];

%T = all_SU{1,1}.xvec;
xvec = all_SU{1,goods(1)}.xvec.CON;
T = find(xvec>=win(1),1):find(xvec<=win(2),1,'last');
ax_time = xvec(T);


Label_vol = {'5-10uL','20uL','40-80uL'};

for a = 1:length(goods)
        clusternames(a) = all_SU{1,goods(a)}.cluster_id;
end


for k = 1:5
    % for i = 1:length(goods)
    %     all_SU{1,goods(i)}.hmat.(string(alignto(k))) = nanzscore(all_SU{1,goods(i)}.hmat.(string(alignto(k))));
    % end

    %% calculate trial-average of z-scored spike rate for each cell
    for i = 1:length(goods)
        this =  all_index{goods(i),4};  %S struct number
        %trials = all_behEvents{1,this}.NoVio;
        %% find trials for each reward vol
        S = all_S{1,this};
        [S.RewardAmount] = convertreward(S.RewardAmount);

        %% find volume response
        %non-violation trials
        
        %remove post-violation trials
        postvio = find(S.vios==1)+1;
        S.vios(postvio) = 1;
        if length(S.vios)>length(S.NoseInCenter)
        S.vios(length(S.NoseInCenter)+1) = [];
        end

        Vol{1,:} = find(S.vios==0 & S.RewardAmount<3 & S.Block==1); %5/10
        Vol{2,:} = find(S.vios==0 & S.RewardAmount==3 & S.Block==1); %2
        Vol{3,:} = find(S.vios==0 & S.RewardAmount>3 & S.Block==1); %40/80

        % Vol{1,:} = find(S.vios==0 & S.RewardAmount<3 ); %5/10
        % Vol{2,:} = find(S.vios==0 & S.RewardAmount==3); %2
        % Vol{3,:} = find(S.vios==0 & S.RewardAmount>3); %40/80

        FVol{1,:} = Vol{1}((Vol{1,:}<length(S.NoseInCenter)))+1;
        FVol{2,:} = Vol{2}((Vol{2,:}<length(S.NoseInCenter)))+1;
        FVol{3,:} = Vol{3}((Vol{3,:}<length(S.NoseInCenter)))+1;

        trials(i,:) = [Vol(1),Vol(2),Vol(3)];
        Ftrials(i,:) = [FVol(1),FVol(2),FVol(3)];

        for v = 1:3
            avgSpk.(alignto{k}){v}(i,:) = mean(all_SU{1,goods(i)}.hmat.(string(alignto(k)))(Vol{v},:), 'omitnan');  % average spike rate for a cell over all trials in a session
            FavgSpk.(alignto{k}){v}(i,:) = mean(all_SU{1,goods(i)}.hmat.(string(alignto(k)))(FVol{v},:), 'omitnan');  % average spike rate for a cell over all trials in a session

        end
    end
end


%% plotting

%sorting
if nargin>5
    if isempty(varargin{1})
    vol5 = avgSpk.(alignto{1}){1, 1};
else
    vol5 = avgSpk.(varargin{1}){1, 1};
    end
else 
    vol5 = avgSpk.(alignto{1}){1, 1};
end

vol5_top = max(vol5,[],2);
for i = 1:length(vol5_top)
    try
    vol5_local(i)= find(vol5(i,:)==vol5_top(i),1);
    catch
        vol5_local(i)=1;
        disp(['has issues sorting', all_index(i,1:3)])
    end
end
vol5_local_sort = [vol5_local;1:length(vol5_top)]';
vol5_local_sort = sortrows(vol5_local_sort);
vol5_indx = vol5_local_sort(:,2);


%Plot
if nargin>6
else
figure('PaperSize',[11 8.5])
for k = 1:5
    for v = 1:3
        p(k) = subplot(3,5,(5*(v-1)+k));
         imagesc(ax_time, linspace(1,size(goods,2)), avgSpk.(string(alignto{k})){v}(vol5_indx,T))

        if v==1
            title(alignto(k))
        end

        clim(caxis)
        xline(0,'--w',LineWidth=1)
        box off
        if k==5
            colorbar
            p(k).Position = [(p(k).Position(1)) (p(k).Position(2)) (p(k-1).Position(3)) (p(k-1).Position(4))];
            drawnow
            yticklabels([])
        elseif k==1
            yticklabels(clusternames(vol5_indx))
            ylabel(Label_vol{v})
        else
            yticklabels([])
        end
    end
    colormap(TriColor)
end

%Population volume encoding on the current trial
figure
for k = 1:5
    subplot(1,5,k)
shadedErrorBar(ax_time,mean(avgSpk.(alignto{k}){1,1}(:,T),1,'omitnan'),sem(avgSpk.(alignto{k}){1,1}(:,T),'omitnan'),'Lineprops',{'Color','b'})
shadedErrorBar(ax_time,mean(avgSpk.(alignto{k}){1,2}(:,T),1,'omitnan'),sem(avgSpk.(alignto{k}){1,2}(:,T),'omitnan'),'Lineprops',{'Color','k'})
shadedErrorBar(ax_time,mean(avgSpk.(alignto{k}){1,3}(:,T),1,'omitnan'),sem(avgSpk.(alignto{k}){1,3}(:,T),'omitnan'),'Lineprops',{'Color','r'})
title(alignto{k})
xlabel("time (s)");xline(0,'--')
box off;set(gca,'TickDir','out');
end
legend('5/10uL', '20uL','40/80uL')
subplot(1,5,1)
ylabel('z-scored ff')
set(gcf,'Color',[1 1 1])
%[h,p]=ttest(avgSpk.Rew{1,1}(:,T)-avgSpk.Rew{1,3}(:,T));%paired
%[h,p]=ttest2(avgSpk.Rew{1,1}(:,T),avgSpk.Rew{1,3}(:,T)); %unpaired
% hold on;
% plot(ax_time(h==1),h(h==1)-0.2,'k*')
sgtitle('Volume encoding on the current trial')

%Population volume encoding on the following trial
figure
for k = 1:5
    subplot(1,5,k)
shadedErrorBar(ax_time,mean(FavgSpk.(alignto{k}){1,1}(:,T),1,'omitnan'),sem(FavgSpk.(alignto{k}){1,1}(:,T),'omitnan'),'Lineprops',{'Color','b'})
shadedErrorBar(ax_time,mean(FavgSpk.(alignto{k}){1,2}(:,T),1,'omitnan'),sem(FavgSpk.(alignto{k}){1,2}(:,T),'omitnan'),'Lineprops',{'Color','k'})
shadedErrorBar(ax_time,mean(FavgSpk.(alignto{k}){1,3}(:,T),1,'omitnan'),sem(FavgSpk.(alignto{k}){1,3}(:,T),'omitnan'),'Lineprops',{'Color','r'})
title(alignto{k})
xlabel("time (s)");xline(0,'--')
box off;set(gca,'TickDir','out');
end
legend('5/10uL', '20uL','40/80uL')
subplot(1,5,1)
ylabel('z-scored ff')
set(gcf,'Color',[1 1 1])
%[h,p]=ttest(avgSpk.Rew{1,1}(:,T)-avgSpk.Rew{1,3}(:,T));%paired
%[h,p]=ttest2(avgSpk.Rew{1,1}(:,T),avgSpk.Rew{1,3}(:,T)); %unpaired
% hold on;
% plot(ax_time(h==1),h(h==1)-0.2,'k*')
sgtitle('Volume encoding on the next trial')


end

% %% mean reward response for each cell
% for a = 1:length(goods)
% figure
% plot(ax_time,avgSpk.Rew{1,1}(a,T),'Color','b'); hold on;
% plot(ax_time,avgSpk.Rew{1,2}(a,T),'Color','k')
% plot(ax_time,avgSpk.Rew{1,3}(a,T),'Color','r')
% ylabel('z-scored ff')
% xlabel("time (s)")
% title('Reward')
% box off
% set(gca,'TickDir','out')
% legend('5/10uL', '20uL','40/80uL')
% end

% % plot the individual neurons
% for i = 1:length(goods)
%     figure('PaperSize',[11 8.5],'Position',[100, 500, 1700,500])
%     for k = 1:5
%         subplot(1,5,k)
%         shadedErrorBar(ax_time,mean(all_SU{1,goods(i)}.hmat.(string(alignto(k)))(trials{i,1},T),'omitnan'),sem(all_SU{1,goods(i)}.hmat.(string(alignto(k)))(trials{i,1},T),'omitnan'),'lineprops',{'color',[0,0,1]})
%         shadedErrorBar(ax_time,mean(all_SU{1,goods(i)}.hmat.(string(alignto(k)))(trials{i,2},T),'omitnan'),sem(all_SU{1,goods(i)}.hmat.(string(alignto(k)))(trials{i,1},T),'omitnan'),'lineprops',{'color',[0,0,0]})
%         shadedErrorBar(ax_time,mean(all_SU{1,goods(i)}.hmat.(string(alignto(k)))(trials{i,3},T),'omitnan'),sem(all_SU{1,goods(i)}.hmat.(string(alignto(k)))(trials{i,1},T),'omitnan'),'lineprops',{'color',[1,0,0]})
%         xline(0,'--');box off; set(gca,'Tickdir','out'); xlabel('time'); ylabel('z-scored ff');axis square; title(alignto{k}); sgtitle(strcat('cluster ',string(all_SU{1,goods(i)}.cluster_id)));
%     end
% end
% 

