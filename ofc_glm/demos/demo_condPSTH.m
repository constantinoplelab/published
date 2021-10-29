%plot conditional PSTHs of neurons

%load data. change datadir on your local machine!
datadir = '~/projects/ofc/data/published/'; %where aggregated data lives
D = load(strcat(datadir,'concatdata_ofc_noOptOut.mat'));
A = D.A;
clear D

set(0,'defaultAxesFontSize',15); %default font size as 15

%%


file_idx = 170; %neuron of choice

%timing vectors
tmesh_start = -2:0.05:4;
tmesh_choice = -4:0.05:4;
tmesh_leavecpoke = -4:0.05:4;

%trial contingincies
left = A{file_idx}.went_right==0;
right = A{file_idx}.went_right==1;
win = A{file_idx}.hits==1;
loss = A{file_idx}.hits==0;
prevwin = A{file_idx}.prev_hits==1;
prevloss = A{file_idx}.prev_hits==0;
V6 = A{file_idx}.chosenval==6 & win;
V12 = A{file_idx}.chosenval==12 & win;
V24 = A{file_idx}.chosenval==24 & win;
V48 = A{file_idx}.chosenval==48 & win;

%combine for easier plotting
volcell = {loss,V6,V12,V24,V48};
volnames = num2str([0,6,12,24,48]');

%% show mean and SEM of each conditional PSTH

%left-right
figure(1)
clf
shadedErrorBar(tmesh_leavecpoke,nanmean(A{file_idx}.hmat_leavecpoke(left,:)), ...
    nanstd(A{file_idx}.hmat_leavecpoke(left,:))/sqrt(sum(left)),...
    'lineprops',{'color','b'});
hold on
shadedErrorBar(tmesh_leavecpoke,nanmean(A{file_idx}.hmat_leavecpoke(right,:)), ...
    nanstd(A{file_idx}.hmat_leavecpoke(right,:))/sqrt(sum(right)),...
    'lineprops',{'color','r'});

vline(0,'k')
xlabel('time to exit center port (s)')
ylabel('rate (Hz)')
title('choice PSTH')
legend('left','right')


%prev.win- prev.loss
figure(2)
clf
shadedErrorBar(tmesh_start,nanmean(A{file_idx}.hmat_start(prevwin,:)), ...
    nanstd(A{file_idx}.hmat_start(prevwin,:))/sqrt(sum(prevwin)),...
    'lineprops',{'color','b'});
hold on
shadedErrorBar(tmesh_start,nanmean(A{file_idx}.hmat_start(prevloss,:)), ...
    nanstd(A{file_idx}.hmat_start(prevloss,:))/sqrt(sum(prevloss)),...
    'lineprops',{'color','r'});
vline(0,'k')
xlabel('time to start (s)')
ylabel('rate (Hz)')
title('reward history PSTH')
legend('prev. win','prev. loss')


%win//loss
figure(3)
clf
shadedErrorBar(tmesh_choice,nanmean(A{file_idx}.hmat_choice(win,:)), ...
    nanstd(A{file_idx}.hmat_choice(win,:))/sqrt(sum(win)),...
    'lineprops',{'color','b'});
hold on
shadedErrorBar(tmesh_choice,nanmean(A{file_idx}.hmat_choice(loss,:)), ...
    nanstd(A{file_idx}.hmat_choice(loss,:))/sqrt(sum(loss)),...
    'lineprops',{'color','r'});
vline(0,'k')
xlabel('time to choice (s)')
ylabel('rate (Hz)')
title('reward PSTH')
legend('win','loss')

%rewarded volume
figure(4)
clf
colormat = linspecer(5,'sequential');
for j = 1:5
    vmask = volcell{j};
    shadedErrorBar(tmesh_choice,nanmean(A{file_idx}.hmat_choice(vmask,:)), ...
        nanstd(A{file_idx}.hmat_choice(vmask,:))/sqrt(sum(vmask)),...
        'lineprops',{'color',colormat(j,:)});
    hold on
end

vline(0,'k')
xlabel('time to choice (s)')
ylabel('rate (Hz)')
title('reward PSTH')
legend(volnames)
