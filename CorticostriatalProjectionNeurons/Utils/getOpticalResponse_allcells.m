function [LR,all_SU,align_spikes]=getOpticalResponse_allcells(all_SU,all_S,all_index, win, which, bins, pt_choice)
% this finds the times of all spikes within the optical stimluation window
% input:
% all_SU = all cells that have optical stimulation
% all_pt = the pulse times of all animals
% all_index = the indexing of SU to pt, since pt is by animal and SU is concatinated
% win = window over which to look for spikes [-0.01 0.1]
% which = 1=1st spike following stimulation
%         2=spikes around the peak response
% bins = time width of the bin defualt 0.001 = 1ms
% pt_choice = which stimulation version do you want
%           12 = 3ms, 200mA, 10Hz stimulation
%           15 = 3ms, 200mA, 20Hz stimulation
%           18 = 3ms, 200mA, 40Hz stimulation

% outputs:
% which = 1
%   LR = struct with:
%     first_latency = latency between light pulse and first spike
%     first_avg_lat = average latency of first spike per cell
%     first_median_lat = median of the latency between the pulse time and first spike per cell
%     first_jitter = std of the latecny to first spike per cell
%     reliability = percentage of pulses that have a spike within a 5ms window around the mean
%     all_spikes = time of all spikes relative to the pulse time
% which = 2
%    LR = struct with; 
%       latency = latency to the peak response to light stimulation after
%       the pulse begins
%       latency_xvec = time window around the light pulse to align hmats to
%       reliability = probability that the neuron will fire around the peak
%       latency
%


% note that this is very similar to getlatency, but has a reduced
% functionality around plotting median and mean options. I just sorta
% settled on mean/std, but median etc could be revisited

%keyboard
% uses 1st 200 pulses, at 3ms, 200mA, 40Hz stimulation
pulsedur = 0.003;
pulseisi = 0.02;
%pt_choice = 18; % 18 for the 40hz, 3ms, 200mA only
%bins = 0.002;   %1ms time bins


if which == 1

    % uses 1st 200 pulses, at 3ms, 200mA, 10Hz stimulation
    pulsenum = 200;
    %win = [0,(pulseisi-pulsedur)];    %window to look at in seconds 0:100ms
    rel_win = 0.005;  % window around the mean in s


    for cells = 1:length(all_SU) % for all cells
        spiketimes = all_SU{cells}.st;
        pt_indx = cell2mat(all_index(cells,4)); %all_indx = ratname, date, clustername, S-index, cluster-index

        pt = all_S{1,pt_indx}.pt;
        try %get pulse times
            pulseON = pt{pt_choice,1}; % 12 for the 10hz, 3ms, 200mA only
        catch
            pulseON = pt{1,pt_choice}; % 12 for the 10hz, 3ms, 200mA only
        end %get pulse times

        for ii = 1:pulsenum %for all pulses
            t1 = pulseON(ii)+pulsedur + win(1,1);    %end of pulse
            t2 = pulseON(ii) + win(1,2);    %start of the next pulse
            these = find(spiketimes>t1 & spiketimes<t2);   %spikes within the time window
            if sum(these)==0        % if there are no spikes
                align_spikes{ii,1} = NaN;
                first_spike(ii) = NaN;    %% not sure if it should be nan or 0.1
            else    % if there are spikes
                align_spikes{ii,1} = spiketimes(these) - pulseON(ii); % zero time to the pulse
                first_spike(ii) = spiketimes(these(1)) - pulseON(ii);  % zero time of first spike to pulse
            end

        end%for all pulses
        all_spikes{cells,1} = cell2mat(align_spikes);   %times of all spikes following the light
        first_latency(:,cells) = first_spike;   %time to the first spike
        first_avg_lat(cells) = mean(first_spike,'omitnan'); %
        first_mid_lat(cells) = median(first_spike,'omitnan');
        first_jitter(cells) = std(first_spike,'omitnan');

        reliability(cells) = length(intersect(find(first_spike>first_avg_lat(cells)-rel_win),find(first_spike<first_avg_lat(cells)+rel_win)))/pulsenum; % percentage of pulses there is a spike  %% fix this

    end % for all cells

    LR = struct('first_latency',first_latency,...
        'first_avg_lat', first_avg_lat,'first_median_lat',first_mid_lat,...
        'first_jitter', first_jitter,'reliability', reliability);

    % plot the latency/reliability/jitter
    figure;
    plot3(LR.first_avg_lat,LR.first_jitter,LR.reliability, '.');
    xlabel('latency')
    ylabel('jitter')
    zlabel('reliability')
    grid on
    %print(gcf,strcat(savedir,'LatJitterRelia.pdf'),'-dpdf','-bestfit');

    %plot the latency
    figure;
    plot(1:length(LR.first_avg_lat),LR.first_avg_lat, '.');
    xlabel('latency')
    ylabel('jitter')
    zlabel('reliability')
    grid on
    % print(gcf,strcat(savedir,'LatRelia.pdf'),'-dpdf','-bestfit');

elseif which == 2
    

    % uses 1st 200 pulses, at 3ms, 200mA, 10Hz stimulation
    pulsenum = 200;
    %win = [0,pulseisi];  %time between the pulse off to the next pulse on

    for cells = 1:length(all_SU) % for all cells
        %time vector
        xvec = win(1,1):bins:win(1,2); 
        d = diff(xvec)/2;
        edges = [xvec(1)-d(1), xvec(1:end-1)+d, xvec(end)+d(end)];



        spiketimes = all_SU{cells}.st;

        pt_indx = all_index{cells,4}; %all_indx = ratname, date, clustername, S-index, cluster-index
        try
        pt = all_S{1,pt_indx}.pt;
        catch
            keyboard
        end
        
        try %get pulse times
            pulseON = pt{pt_choice,1}; % for the 10hz, 3ms, 200mA only start of pulse
        catch
            pulseON = pt{1,pt_choice}; % for the 10hz, 3ms, 200mA only
        end %get pulse times


        for ii = 1:pulsenum %for all pulses
            t1 = pulseON(ii) + win(1,1);    %start of the pulse + window
            t2 = pulseON(ii) + win(1,2);    %start of the next pulse
            these = find(spiketimes>t1 & spiketimes<t2);   %spikes within the time window
            %zero the spike times to the end of the pulse
            if sum(these)==0        % if there are no spikes
                align_spikes{ii,cells} = NaN;
                spike_times(these) = NaN;
                hmat(ii,1:length(xvec)) = 0;  %% maybe this should be nan?

            else    % if there are spikes
                align_spikes{ii,cells} = spiketimes(these) - pulseON(ii); % zero time to the start of the pulse
                spike_times(these) = spiketimes(these) - pulseON(ii);
 

                [n,edges] = histcounts(spike_times(these), edges);
                hmat(ii,:) = n; %number of spikes n/.b = frequency at that bin

            end

        end%for all pulses
        all_SU{1,cells}.phmat = hmat;    %save the hmat of the pulse responses
        all_SU{1,cells}.latency_xvec = xvec;
        psth(cells,:) = mean(hmat);

        zero_location = find(xvec==0);         %find where zero is
        time2peak(cells) = find(psth(cells,zero_location:end)==max(psth(cells,zero_location:end)),1)+zero_location-1; %time bin with the highest density of spikes
        % reliability(cells) = psth(cells,time2peak(cells));    %mean number of spikes at time2peak

        reliability(cells)=  sum(hmat(:,time2peak(cells))>0)/pulsenum;    %number of pulses with a spike at latency/#of pulses


    end % for all cells

    latency(1,:) = xvec(time2peak);    %bin of peak time x bin duration x 1000ms/s

    LR = struct('latency',latency,'latency_indx',time2peak,'xvec',xvec,'reliability', reliability);

    %plot the latencyxreliability
    figure;
    plot(LR.latency*1000, LR.reliability*100, '.');
    xlabel('latency (ms)')
    ylabel('reliability (%)')
    xlim([-2 22]);
    ylim([-10 100])
    yticks(0:20:100)
    grid off
    box off
    ax = gca;
    ax.TickDir = 'out';
    ax.Color = 'none';
    if pt_choice ==12
        title('10Hz')
    elseif pt_choice ==15
        title('20Hz')
    elseif pt_choice == 18
        title('40Hz')
    end

else
    keyboard %you put in a which value that does not exist
end


