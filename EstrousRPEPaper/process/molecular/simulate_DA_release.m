function simulate_DA_release
%Dopamine simulation for Golden et al. (Figure 5g-h). 

T = -.05:.001:.100; %time vector
tstart = find(T==0);
rew = 1:3;
mycolors = [0 0 1; .55 0 .8; 1 0 0];
gammas = [.7 .85]; %terms specifying reuptake
timeStepS = 0.001;
durationS = .015;
ntrials = 500;
evoked_FRs = [0 1.5 3]; %evoked DA firing rate for different RPEs.
spont = 15; %spontaneous firing rate.
figure
ctr = 1;

z = zeros(ntrials, length(T), length(rew));

for k = 1:length(rew)
    
    epsilon = randn(ntrials,length(T))./100;
    
    vt = rand(size(T))./20;
    
    spikesPerS = evoked_FRs(k)*20;
    
    
    [spontspikes, ~] = makeSpikes(timeStepS, spont, length(T), ntrials);
    z(:,:,k) = spontspikes(:,1:length(T));
    
    [spikes, ~] = makeSpikes(timeStepS, spikesPerS, durationS, ntrials);
    
    
    if k==1
        z(:,tstart:tstart+length(spikes(1,:))-1,k) = zeros(ntrials,length(tstart:tstart+length(spikes(1,:))-1));
    else
        z(:,tstart:tstart+length(spikes(1,:))-1,k) = z(:,tstart:tstart+length(spikes(1,:))-1, k)+spikes;
    end
end

y = nan(ntrials, length(T), 1);
c = y;



examp = nan(2, length(T));
for M =1:2 %loop over gammas
    Y = [];

    %loop over rewards
    for k = 1:length(rew)
        y(:,1) = 0;
        c(:,1) = 0;
        for ll = 1:ntrials
            for t = 2:length(T)
                %convolve with each trial
                c(ll,t) = gammas(M)*c(ll,t-1) + z(ll,t,k);
                y(ll,t) = c(ll,t) + epsilon(ll,t);
            end
        end
        
      Y = [Y; mean(y)]; 
      yy = y';
      yy = yy(:);
    end
    
    ix2 = find(sum(spikes')==3);
    examp(M,:) =  y(ix2(1),:);
    
    subplot(2,3,ctr)
    plot(T, examp(M,:), 'Color', [0 0 0]); hold on
    set(gca, 'TickDir', 'out'); box off
    yl = ylim;
    xlim([-.02 .05]);
    ylim([-.1 2]);
    ctr = ctr+1;
    
    subplot(2,3,ctr)
    for j = 1:length(Y(:,1))
        plot(T, Y(j,:), 'Color', mycolors(j,:)); hold on
    end
    set(gca, 'TickDir', 'out'); box off
    ylim([-.05 .5]);
    xlim([-.02 .05]);
    line([0 0], [-.1 .5], 'Color', [0 0 0], 'LineStyle', '--');
    ylabel('normalized [DA]');
    xlabel('Time from phasic spikes');
    title(strcat(['Gamma = ', num2str(gammas(M))]));
    ctr = ctr+1;
    set(gcf, 'Color', [1 1 1]);

    subplot(2,3,3);
    mc = [0 0 0; .5 .5 .5];
    auc = nan(length(rew),1);
    for k = 1:length(rew)
        t1 = find(T==0);
        t2 = find(abs(T-0.02)<eps);
        baseline = mean(Y(k,1:t1-1)');
        Y(k,:) = Y(k,:)-baseline;
        auc(k) = trapz(t1:t2,Y(k,t1:t2));
    end
    plot(1:length(rew), auc, 'Color', mc(M,:)); hold on
    plot(1:length(rew), auc, 'Color', mc(M,:), 'Marker', '.', 'MarkerSize', 10); hold on
    xlabel('Reward');
    ylabel('AUC of phasic DA');
    set(gca, 'TickDir', 'out'); box off;
    xlim([.5 length(rew)+0.5]);
    ctr = ctr+1;
end

end

function [spikes, ix] = makeSpikes(timeStepS, spikesPerS, durationS, numTrains)

if (nargin < 4)
    numTrains = 1;
end
times = [0:timeStepS:durationS];
spikes = zeros(numTrains, length(times));
for train = 1:numTrains
    vt = rand(size(times));
    spikes(train, :) = (spikesPerS*timeStepS) > vt;
    ix{train} = find(spikes(train,:)==1);
end
end