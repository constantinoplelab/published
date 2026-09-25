function [highc,highe,lowc,lowe,mixc,mixe]= wt_notsided(ExpA,ConA,Raw_data,varargin)
%no side bias wt
% Inputs: 
    % ExpA: Opto A struct
    % ConA: control A struct
    % Rawdata: Both control and opto S structs
    % varargin: 1 - plots each animal wt curves

for a = 1:length(ExpA)
    %generate the mean wt across all sessions
    rat = ConA{a}.ratname;
    S.pd = [Raw_data.(rat).exp.pd,Raw_data.(rat).con.pd];
    S.peh = [Raw_data.(rat).exp.peh,Raw_data.(rat).con.peh];
    A = parse_data_from_mysql(S,[],1);

    %threshdold all waittimes
    AA = load(fullfile('Z:\ProcessedRatData\A_Structs',['ratTrial_',ConA{a}.RatName,'.mat']));
    wait_thresh = AA.A.wait_thresh;
    disp([ConA{a}.RatName, string(wait_thresh)])
    A.wait_time(A.wait_time>wait_thresh) = NaN;
    ConA{a}.wait_time(ConA{a}.wait_time>wait_thresh) = NaN;
    ExpA{a}.wait_time(ExpA{a}.wait_time>wait_thresh) = NaN;

    %detrend all sessions and individual sessions
    A = detrendwt(A);
    ConA{a} = detrendwt(ConA{a});
    ExpA{a} = detrendwt(ExpA{a});

    [hi, low, mix] = blocks(A);
    
    %normalize the wt by the mean mix 20uL wt
    [hic, loc, mc] = blocks(ConA{a});
    [hie, loe, me] = blocks(ExpA{a});

    % highc(a,:) = hic.wt/mix.wt(3);
    % highe(a,:) = hie.wt/mix.wt(3);
    % 
    % lowc(a,:) = loc.wt/mix.wt(3);
    % lowe(a,:) = loe.wt/mix.wt(3);
    % 
    % mixc(a,:) = mc.wt/mix.wt(3);
    % mixe(a,:) = me.wt/mix.wt(3);


    % normalize by control mix
    highc(a,:) = hic.wt/mc.wt(3);
    highe(a,:) = hie.wt/mc.wt(3);

    lowc(a,:) = loc.wt/mc.wt(3);
    lowe(a,:) = loe.wt/mc.wt(3);

    mixc(a,:) = mc.wt/mc.wt(3);
    mixe(a,:) = me.wt/mc.wt(3);


% if ~isempty(varargin)
%     figure
%     subplot(1,2,1)
%     plot(1:5,highc(a,:),'r');hold on
%     plot(1:5,highe(a,:),'color',[1,0.5,0.5])
%     plot(1:5,lowe(a,:),'color',[0.5,0.5,1])
%     plot(1:5,lowc(a,:),'color',[0,0,1])
%     box off; set(gca,'TickDir','out');xlabel('Offered Reward');xticks(1:5); ylabel('normalized WT');
%     ylim([0.8 1.3])
%     subplot(1,2,2)
%     plot(1:5,mixc(a,:),'color',[0,0,0]);hold on
%     plot(1:5,mixe(a,:),'color',[0.5,0.5,0.5])
%      box off; set(gca,'TickDir','out');xlabel('Offered Reward');xticks(1:5);
%       ylim([0.8 1.3])
%      sgtitle(rat)
%      legend('control','opto',Location='southeast')
% end

if ~isempty(varargin)
    figure
    subplot(1,2,1)
    % plot(1:5,highc(a,:),'r');hold on
    %     plot(1:5,mixc(a,:),'color',[0,0,0]);hold on
    % plot(1:5,lowc(a,:),'color',[0,0,1])
    plot(1:5,highc(a,:),'r');hold on
    plot(1:5,highe(a,:),'color',[1,0.5,0.5]); hold on
    plot(1:5,lowc(a,:),'color',[0,0,1])
    plot(1:5,lowe(a,:),'color',[0.5,0.5,1])
    
    box off; set(gca,'TickDir','out');xlabel('Offered Reward');xticks(1:5); ylabel('normalized WT');
    ylim([0.8 1.3])
    subplot(1,2,2)
    % plot(1:5,highe(a,:),'color',[1,0.5,0.5]); hold on
    % plot(1:5,mixe(a,:),'color',[0.5,0.5,0.5])
    %     plot(1:5,lowe(a,:),'color',[0.5,0.5,1])
    plot(1:5,mixc(a,:),'color',[0,0,0]);hold on
    % plot(1:5,lowc(a,:),'color',[0,0,1])
    plot(1:5,mixe(a,:),'color',[0.5,0.5,0.5])
    % plot(1:5,lowe(a,:),'color',[0.5,0.5,1])

     box off; set(gca,'TickDir','out');xlabel('Offered Reward');xticks(1:5);
      ylim([0.8 1.3])

     sgtitle(rat)
     legend('control','opto',Location='southeast')
end

end

if ~isempty(varargin)
figure
subplot(1,3,1)
shadedErrorBar(1:5, mean(mixc,'omitnan'),sem(mixc,'omitnan'),'lineprops',{'color',[0,0,0]})
shadedErrorBar(1:5, mean(mixe,'omitnan'),sem(mixe,'omitnan'),'lineprops',{'color',[0.5,0.5,0.5]})
hold on;
h = ttest(mixe-mixc);
%normal =  kstest(mixe(:,4)-mixc(:,4));
plot(find(h==1),h(h==1)+0.3,'*k')

ylabel('Nomarlized wait time')
xticks(1:5); xticklabels([5,10, 20, 40, 80]); xlim([0 6]); xlabel('Reward Volume (uL)')
set(gca,'TickDir','out');box off;
% ylim([0.88 1.32]),yticks(0.9:0.1:1.3)

subplot(1,3,2)
shadedErrorBar(1:5, mean(highc,'omitnan'),sem(highc,'omitnan'),'lineprops',{'color',[1,0,0]})
shadedErrorBar(1:5, mean(highe,'omitnan'),sem(highe,'omitnan'),'lineprops',{'color',[1,0.5,0.5]})
hold on;

% change this to a permutation test?
% figure
% histogram(highe);
% normal = kstest(highe-highc);
% Last tests showed these were normally distributed.

[h,p] = ttest(highe-highc);
disp(['high',string(p)])
p = permute_test(highe,highc,100);
disp(['high permute',string(p)])
plot(find(h==1),h(h==1)+0.3,'*k')

shadedErrorBar(1:5, mean(lowc,'omitnan'),sem(lowc,'omitnan'),'lineprops',{'color',[0,0,1]})
shadedErrorBar(1:5, mean(lowe,'omitnan'),sem(lowe,'omitnan'),'lineprops',{'color',[0.5,0.5,1]})
hold on;
[h,p] = ttest(lowe-lowc);
disp(['low',string(p)])
p = permute_test(lowe,lowc,100);
disp(['high permute',string(p)])
plot(find(h==1),h(h==1)+0.3,'*k')

ylabel('Nomarlized wait time')
xticks(1:5); xticklabels([5,10, 20, 40, 80]); xlim([0 6]); xlabel('Reward Volume (uL)')
set(gca,'TickDir','out');box off;
ylim([0.88 1.32]),yticks(0.9:0.1:1.3)


subplot(1,3,3)

lowchange = mean(lowe-lowc,2,'omitnan');
highchange = mean(highe-highc,2,'omitnan');
mixchange = mean(mixe-mixc,2,'omitnan');

plot(1:3,[ lowchange,mixchange,highchange]','color',[0.5 0.5 0.5])
hold on;
plot(2,mixchange,'.','color',[0,0,0]);
plot(2.25, mean(mixchange,'omitnan'),'_k')
errorbar(2.25,mean(mixchange,'omitnan'),sem(mixchange,'omitnan'),'-k')

plot(1,lowchange,'.','color',[0,0,1]);
plot(1.25, mean(lowchange,'omitnan'),'-b')
errorbar(1.25,mean(lowchange,'omitnan'),sem(lowchange,'omitnan'),'_b')

plot(3,highchange,'.','color',[1,0,0]);
plot(3.25, mean(highchange,'omitnan'),'-r')
errorbar(3.25,mean(highchange,'omitnan'),sem(highchange,'omitnan'),'_r')
box off
xticklabels({'low','mix','high'}), set(gca,'TickDir','out')
yline(0,'--');xlim([0.5, 3.5]); ylabel('Delta WT by opto'); xlabel('Block');

% keyboard
% 
% %plot(1:3,[mixe(:,3)-mixc(:,3),lowe(:,3)-lowc(:,3),highe(:,3)-highc(:,3)],"Color",[0.5,0.5,0.5]); hold on;
% 
% %plot(1,mixe(:,3)-mixc(:,3),'.','color',[0,0,0]);
% plot(1.25, mean(mixe(:,3)-mixc(:,3),'omitnan'),'_k')
% errorbar(1.25,mean(mixe(:,3)-mixc(:,3),'omitnan'),sem(mixe(:,3)-mixc(:,3),'omitnan'),'-k')
% 
% plot(2,lowe(:,3)-lowc(:,3),'.','color',[0,0,1]);
% errorbar(2.25,mean(lowe(:,3)-lowc(:,3),'omitnan'),sem(lowe(:,3)-lowc(:,3),'omitnan'),'-b')
% plot(2.25,mean(lowe(:,3)-lowc(:,3),'omitnan'),'_b')
% 
% plot(3,highe(:,3)-highc(:,3),'.','color',[1,0,0]);
% errorbar(3.25,mean(highe(:,3)-highc(:,3),'omitnan'),sem(highe(:,3)-highc(:,3),'omitnan'),'-r')
% plot(3.25, mean(highe(:,3)-highc(:,3),'omitnan'),'_r');
% 
% xlim([0 4]); xticks(1:4); xticklabels({'mix', 'low','high'}); xlabel('block'); 
% ylabel('Change in 20uL wait time'); yline(0,'--')
% set(gca,'TickDir','out');box off;
% 
% [h,p] = ttest2(lowe(:,3)-lowc(:,3),highe(:,3)-highc(:,3));
end