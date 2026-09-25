function [dprime,sig,dvar] = dprimespecifyvol(goods,all_SU,all_S,all_index,win,ind,V1,V2,plotit,event)
%% d' of 40/80
datpath = 'C:\Users\mld9131\Documents\GitHub\Papers\MLD\Figure2\Data';
load(fullfile(datpath,'TriColor.mat'))

for i = 1:length(goods)
S = all_S{all_index{goods(i),4}};
S.RewardAmount = convertreward(S.RewardAmount);
clusternames(i) = all_SU{goods(i)}.cluster_id;

%find volume trials
Vol{1} = find(S.RewardAmount==V1 & S.vios==0 & S.Block==1);
Vol{2} = find(S.RewardAmount==V2 & S.vios==0 & S.Block==1);
% Vol{1} = find(S.RewardAmount==V1 & S.vios==0);
% Vol{2} = find(S.RewardAmount==V2 & S.vios==0);
trials(i,:) = [Vol(1),Vol(2)];
end %cells

[dprime,sig] = vol_dprime(trials,goods,all_S,all_SU,all_index,event,[-1 3],'signed',200);

signed = 'signed';
[d,dvar,dprime2,~] = makeDprime(trials, goods, all_S,all_SU,all_index,event,0);

% %Convert NANs to zero
c = isnan(dprime2);
dprime2(c) = 0;

if plotit
xvec = all_SU{goods(1)}.xvec.CON;
T = find(xvec>=win(1),1):find(xvec<=win(2),1,'last');
figure
tiledlayout(1,2)
nexttile
imagesc(xvec(T), 1:length(goods),dprime2(ind,T))
yticks(1:length(goods));yticklabels(goods(ind))
colorbar;colormap(TriColor);clim([-2 2])
% title('dprime 5/10') 
title(strcat('dprime ',string(V1),'/', string(V2)))
set(gca,'TickDir','out','box','off')
xline(0,'--');xlabel('Time (s)'); xticks(-1:3)

nexttile
imagesc(xvec(T), 1:length(goods),d(ind,T))
yticks(1:length(goods));yticklabels(goods(ind))
colorbar;colormap(TriColor);clim([-10 10])
% title('FR difference 5/10') 
title(strcat('FR difference ',string(V1),'/',string(V2)))
xline(0,'--');xlabel('Time (s)'); xticks(-1:3)

end






