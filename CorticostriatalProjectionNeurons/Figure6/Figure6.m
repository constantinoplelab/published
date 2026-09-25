function Figure6(datpath)

%% Figure 6: 
% Dmixed PCA shows that OFC-CPi perturbation disrupts block encoding 
% but not reward encoding. 

%% A & B Schematics
%% C Top Block dPC
plotFunction = @dpca_plot_default;
figure('color','w','Units','centimeters','Position',[1 1 18 28])
tiledlayout(6,3)

event = 'Rew';
figurename = 'L076_29-Jan-2026';
% datpath = 'C:\Users\mld9131\Documents\GitHub\Papers\MLD\Figure6\Data';
load(fullfile(datpath,[figurename,'.mat']))
c = 'High-mix';
% Block PCs

    %Control Top Block dPC high-mix
nexttile
        [Zfull_control] = prepplot(FR,event,W);
        vcolors = [0 .4 0; 1 .65 0; 1 0 0]; %green, orange, red control
        PC = 11;
        ywin = [-50 50];
        plotFunction(Zfull_control(PC,:,:,:,:), time, ywin, ...
                [], 4, 0, ...
                [], [],vcolors)
        xline(0,'--')
        title('Top block dPC','Control')
        ylim([-25 25])
        ylabel('')
         legend('20uL mix','20uL high','40uL mix','40uL high','80uL mix','80uL high','Location','eastoutside')

    %Opto Top Block dPC
nexttile        
        [Zfull_opto] = prepplot(FRo,event,W);
         vcolors =  [0.4 .8 0.4; 1 .8 .5; 1 0.5 0.5]; % green. orange, red opto
        ywin = [-50 50];
        plotFunction(Zfull_opto(PC,:,:,:,:), time, ywin, ...
                [], 4, 0, ...
                [], [],vcolors)
        xline(0,'--')
        title('Top block dPC','Opto')
        ylim([-25 25])


    %Session Average Difference
nexttile
% load(fullfile('\\constantinoplelab.cns.nyu.edu\server2\PhysiologyData\Maggie\Chronic_implant\Npxl\Optotagged_cells\figures\dPCA\',c,'\Data-rew\Data','AverageSessionddPCA.mat'))
load(fullfile(datpath,'High-mix_AverageSessionddPCA.mat'))    
        shadedErrorBar(time, mean(D_PSTH_block_avg), sem(D_PSTH_block_avg))
        shadedErrorBar(time, mean(D_PSTH_block_opto_avg), sem(D_PSTH_block_opto_avg),'lineprops',{'color',[0.5 0.5 0.5]});
        hold on
        h = ones(1,length(time));
        plot(time(sig_block),h(sig_block)+19,'.k')
        title('Block difference across sessions')
        ylabel([ '\Delta' 'FR'])
        ylim([-5 20])
        xlim([-1 3])

%% D Reward PCs
    %Control Top Reward dPC
nexttile

        [Zfull_control] = prepplot(FR,event,W);
        vcolors = [0 .4 0; 1 .65 0; 1 0 0]; %green, orange, red controlPC = 3;
        ywin = [-50 50];
        PC = 3;
        plotFunction(Zfull_control(PC,:,:,:,:), time, ywin, ...
                [], 4, 0, ...
                [], [],vcolors)
        xline(0,'--')
        title('Top Reward dPC','Control')
            ylim([-75 100])
        ylabel('')
         legend('20uL mix','20uL high','40uL mix','40uL high','80uL mix','80uL high','Location','eastoutside')

        
    % Opto Top Reward dPC
nexttile

        [Zfull_control] = prepplot(FRo,event,W);
        
         vcolors =  [0.4 .8 0.4; 1 .8 .5; 1 0.5 0.5]; % green. orange, red opto
        ywin = [-50 50];
        plotFunction(Zfull_control(PC,:,:,:,:), time, ywin, ...
                [], 4, 0, ...
                [], [],vcolors)
        xline(0,'--')
        title('Top Reward dPC','Opto')
            ylim([-75 100])

    %Session Average Difference
nexttile    
        shadedErrorBar(time, mean(D_PSTH_Rew_avg), sem(D_PSTH_Rew_avg))
        shadedErrorBar(time, mean(D_PSTH_Rew_opto_avg), sem(D_PSTH_Rew_opto_avg),'lineprops',{'color',[0.5 0.5 0.5]});
        hold on
        h = ones(1,length(time));
        plot(time(sig_rew),h(sig_rew)+68,'.k')
        title('Reward difference across sessions')
        ylabel([ '\Delta' 'FR'])
                   ylim([-10 70])
            xlim([-1 3])


for s = 1:6
    nexttile(s)
    xlim([-1 3])
    set(gca, 'TickDir','out','Box','off')
    xline(0,'--')
end


%% E: Low Mix Block
%% Low-mix
% load 1/13/26 dpca
figurename = 'L076_13-Jan-2026';
% load(fullfile('\\constantinoplelab.cns.nyu.edu\server2\Physiologydata\Maggie\Chronic_implant\Npxl\Optotagged_cells\figures\dPCA\Low-mix\Data-rew\Data',[figurename,'.mat']))
load(fullfile(datpath,[figurename,'.mat']))
c = 'Low-mix';

% Block PCs
    %Control Top Block dPC low-mix
nexttile
[Zfull_control] = prepplot(FR,event,W);
        vcolors = [.6 .2 .9; .1 .5 .9; 0 .4 0]; %purple, blue, green control
        PC = 10;
        ywin = [-50 50];
        plotFunction(Zfull_control(PC,:,:,:,:), time, ywin, ...
                [], 4, 0, ...
                [], [],vcolors)
        xline(0,'--')
        title('Top block dPC','Control')
        ylim([-25 25])
        ylabel('')
        legend('5uL mix','5uL low','10uL mix','10uL low','20uL mix','20uL low','Location','eastoutside')


    % Opto Top Block dPC
nexttile
        [Zfull_opto] = prepplot(FRo,event,W);
        vcolors = [.8 0.6 1;  .3 .7 .9; 0.25 .7 0.25]; % purple, blue, green opto
        ywin = [-50 50];
        plotFunction(Zfull_opto(PC,:,:,:,:), time, ywin, ...
            [], 4, 0, ...
            [], [],vcolors)
        xline(0,'--')
        title('Top block dPC','Opto')
        ylim([-25 25])

    %Session Average Difference
nexttile

        % load(fullfile('\\constantinoplelab.cns.nyu.edu\server2\PhysiologyData\Maggie\Chronic_implant\Npxl\Optotagged_cells\figures\dPCA\',c,'\Data-rew\Data','AverageSessionddPCA.mat'))
    load(fullfile(datpath,'Low-mix_AverageSessionddPCA.mat'))
        shadedErrorBar(time, mean(D_PSTH_block_avg), sem(D_PSTH_block_avg))
        shadedErrorBar(time, mean(D_PSTH_block_opto_avg), sem(D_PSTH_block_opto_avg),'lineprops',{'color',[0.5 0.5 0.5]});
        hold on
        h = ones(1,length(time));
        plot(time(sig_block),h(sig_block)+19,'.k')
        title('Block difference across sessions')
        ylabel([ '\Delta' 'FR'])
        ylim([-5 20])
        xlim([-1 3])


%% F. Reward PCs
        %Control Top Reward dPC
nexttile
            PC = 3;
            vcolors = [.6 .2 .9; .1 .5 .9; 0 .4 0]; %purple, blue, green control
            plotFunction(Zfull_control(PC,:,:,:,:), time, ywin, ...
                [], 4, 0, ...
                [], [],vcolors)
            xline(0,'--')
            title('Top Reward dPC', 'Control')
            ylim([-75 100])
            ylabel('')
            legend('5uL mix','5uL low','10uL mix','10uL low','20uL mix','20uL low','Location','eastoutside')

        % Opto Top Reward dPC
nexttile
            vcolors = [.8 0.6 1;  .3 .7 .9; 0.25 .7 0.25]; % purple, blue, green opto
            plotFunction(Zfull_opto(PC,:,:,:,:), time, ywin, ...
                [], 4, 0, ...
                [], [],vcolors)
            xline(0,'--')
            xlim([-1 3])
            title('Top Reward dPC','Opto')
            ylim([-75 100])

        %Session Average Difference
nexttile          
            shadedErrorBar(time, mean(D_PSTH_Rew_avg), sem(D_PSTH_Rew_avg))
            shadedErrorBar(time, mean(D_PSTH_Rew_opto_avg), sem(D_PSTH_Rew_opto_avg),'lineprops',{'color',[0.5 0.5 0.5]});
            hold on
            h = ones(1,length(time));
            plot(time(sig_rew),h(sig_rew)+62,'.k')
            title('Reward difference across sessions')
            ylabel([ '\Delta' 'FR'])

            ylim([-10 70])
            xlim([-1 3])

%% G.
open(fullfile(datpath,'Figure6\Data\delta_FR_during_opto.fig"'))
