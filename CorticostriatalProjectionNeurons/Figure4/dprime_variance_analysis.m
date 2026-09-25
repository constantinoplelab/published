load(fullfile('\\constantinoplelab.cns.nyu.edu\server2\PhysiologyData\Maggie\Chronic_implant\Npxl\Optotagged_cells','OFC-DS_projection_neurons_non-Stimulated.mat'));
DLS = 1:length(all_SU_DS);

xvec_200 = -4:0.2:8;
xvec_100 = -4:0.1:8;
xvec = -4:0.05:8;
%get d' of DS projectors usign 200ms bins
win = [-2,6];
trials_DS = volbin_all_cells(all_SU_DLS_200,DLS,all_S_DS,all_index_DS,win,'Rew',0);

signedDprime = 'signed';
event = 'Rew';
% event = 'Son';

%% Compare variance of reward encoding
[dprime_DS_new,sig_DS_new,~,dvar] = vol_dprime(trials_DS(:,[1,3]),DLS,all_S_DS,all_SU_DLS_200,all_index_DS,event,win,signedDprime,shufnum); %dprime
T = find(xvec_200>=win(1),1):find(xvec_200<=win(2),1,'last');
times = xvec_200(T);
c = isnan(dprime_DS);
d = dprime_DS;
d(c) = 0;
ind_DS = plotdprime(d,T,xvec_200,0,1);
DLS_sort = DLS(ind_DS);


for b = 1:5
z = 5;  %number of groups
%b = 1;  %which group to study
a = (length(DLS_sort)/z);
section  = round(1+a*(b-1):a*(b-1)+a); 
value = DLS_sort(section);

win = [0 1];
T = find(xvec_200>=win(1),1):find(xvec_200<=win(2),1,'last');
dmeanvar(:,b) = mean(dvar(value,T),2,'omitnan'); %variance in encoding of 5/10 & 40/80
end

Anv = anova(dmeanvar);
% Anv = anova1(dmeanvar); % does the same thing as above
Stat = stats(Anv);
figure
% plot(g1+0.1, y,'.k')
hold on
errorbar(mean(dmeanvar),sem(dmeanvar),'_k')
set(gca,'box','off','TickDir','out')
xticks(1:5); xlim([0.5,5.5]); xlabel('Quintile (Postive encoders -> Negative encoders)'); 
ylabel("Mean variance of 5/10 & 40/80 encoding")
set(gcf,'color','w')
title(strcat('1-way ANOVA, p=',string(Stat.pValue(1))))
exportgraphics(gca,fullfile('Z:\Maggie\Papers\Physiology\Data\Rebuttle','pooled_variance_quintiles.pdf'))



%% Compare variance of reward encoding
[dprime_DS_new,sig_DS_new,~,dvar] = vol_dprime(trials_DS(:,[1,3]),DLS,all_S_DS,all_SU_DLS_200,all_index_DS,event,win,signedDprime,shufnum); %dprime
DLS_sort = DLS(ind_DS);


for b = 1:5
z = 5;  %number of groups
%b = 1;  %which group to study
a = (length(DLS_sort)/z);
section  = round(1+a*(b-1):a*(b-1)+a); 
value = DLS_sort(section);


win = [-2,6];
T = find(xvec_100>=win(1),1):find(xvec_100<=win(2),1,'last');
[pvalued4080,psig4080,dvar4080] = dprimespecifyvol(value,all_SU_DLS_100,all_S_DS,all_index_DS,[-1 3],1:length(value),4,5,1);
[pvalue1020,psig1020,dvar510] = dprimespecifyvol(value,all_SU_DLS_100,all_S_DS,all_index_DS,[-1 3],1:length(value),1,2,1);

win = [0 1];
T = find(xvec_100>=win(1),1):find(xvec_100<=win(2),1,'last');
dmeanvar4080(:,b) = mean(dvar4080(:,T),2,'omitnan'); %variance in encoding of 5/10 & 40/80
dmeanvar510(:,b) = mean(dvar510(:,T),2,'omitnan'); %variance in encoding of 5/10 & 40/80
end

g1 = [ones(12,1);ones(12,1)*2;ones(12,1)*3;ones(12,1)*4;ones(12,1)*5;ones(12,1);ones(12,1)*2;ones(12,1)*3;ones(12,1)*4;ones(12,1)*5];
g2 = [ones(60,1); ones(60,1)*2]; %1 = 5/10, 2 = 40,80
meanvar510 = reshape(dmeanvar510,60,1);
meanvar4080 = reshape(dmeanvar4080,60,1);
y = [meanvar510;meanvar4080];
tbl = table(g1,g2,y);
Anv = anova(tbl,'y ~ g1 + g2+ g1:g2'); %two-way anova
Stat = stats(Anv);


figure('color','w')
% plot(g1+0.1, y,'.k')
% hold on
errorbar(mean(dmeanvar510),sem(dmeanvar510),'_b')
hold on
errorbar((1:5)+0.1,mean(dmeanvar4080),sem(dmeanvar4080),'_r')

set(gca,'box','off','TickDir','out')
xticks(1:5); xlim([0.5,5.5]); xlabel('Quintile (Postive encoders -> Negative encoders)'); 
ylabel("Mean variance")
set(gcf,'color','w')
title(strcat('2-way ANOVA, p=',string(Stat.pValue(3))))
legend('510','4080')
exportgraphics(gca,fullfile('Z:\Maggie\Papers\Physiology\Data\Rebuttle','individual_variance_quintile.pdf'))
