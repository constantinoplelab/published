function [] = graphgen_PAIRS(pairs,pairs_ref,meanvals_dat,meanvals)
%plots angle distributions and visualizes PAIRS statistic
%
%
%written by dlh Oct. 2021


figure(8)
clf
histogram(pairs_ref,100,'Normalization','probability')
CI_95 = 1.95*std(pairs_ref);
vline(mean(pairs_ref)+CI_95,'k')
vline(mean(pairs_ref)-CI_95,'k')
vline(pairs,'r')
xlabel('PAIRS')
title('Reference data PAIRS comparison to data PAIRS statistic (red)')
set(gca,'fontsize',15)

%% make a simplified plot showing angle distribution for real data vs grand reference data

figure(3)
clf
hold on

bins = linspace(0,pi/3,100);
histogram(meanvals(:),bins,'Normalization','probability')
hold on
histogram(meanvals_dat,bins,'Normalization','probability')
title('Distribution of angles for population, PSTH clustering')
xlabel('angle (rad)')
xlim([0,pi/3])
ylabel('p')
legend('reference data','true data')

vline(median(meanvals_dat),'k')
vline(median(meanvals(:)),'k')
set(gca,'fontsize',15)


