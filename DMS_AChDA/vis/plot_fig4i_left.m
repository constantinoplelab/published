function plot_fig4i_left(datadir)
% Plot reward port poke probability during the delay period on trials
% following hit trials with short vs long delay. Mixed block only.
% Short and long trials correspond to bottom and top quartiles of
% reward delay distribution. N = 16 rats.
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.

load(fullfile(datadir, 'data-published/ratlist.mat'), 'ratList');

conditions = {'post short delay', 'post long delay'};
pokeProb = cell(length(ratList), length(conditions));

for rr=1:length(ratList)
    rat = ratList{rr};
    fprintf('%i of %i: %s\n', rr, length(ratList), rat)
    
    [pokeProb(rr,:),bins] = pokePbyPrevDelay(datadir, rat, conditions);
   
end

pokeProb_ = cell(1,length(conditions));
for cc=1:length(conditions)
    pokeProb_{cc} = cell2mat(pokeProb(:,cc));
end


colors = getColorScheme('delay');
mycolor{1} = colors{1};
mycolor{2} = colors{end};

figure;
set(gcf, units='inches', position=[9,4,2.2,1.5])
for cc=1:length(conditions)
    plotPretty(log2(bins), pokeProb_{cc}, mycolor{cc})
    hold on
end
xlabel('Time in trial (s, log_{2})')
xlim([1 4])
ylabel('P(rew port poke)')
set(gca, box='off', tickdir='out', fontsize=8)
legend(conditions, box='off')

end

function [pokeProb_, bins] = pokePbyPrevDelay(datadir, ratname, conditions)

file = strcat('SStruct_', ratname, '.mat');
load(fullfile(datadir, 'data-published\S_Structs', file), 'S');
S_final = get_finalTrainingStageSessions(S);

dt = 0.3;
maxT = 100.; 
bins = 0.5:dt:maxT;

N = zeros(size(bins)); % total number of pokes per time bin
delay_q = get_delayQuartile(datadir, ratname);
lb = delay_q(2);
ub = delay_q(end-1);

for sess=1:length(S_final.pd)
    spd = S_final.pd{sess};
    spd.RewardDelay(spd.RewardDelay<0.75) = nan;
    speh = S_final.peh{sess};
    posthit = [0; spd.hits(1:end-1)];
    postvio = [0; spd.vios(1:end-1)];
    rd_prev = [nan; spd.RewardDelay(1:end-1)];
    for cond=1:length(conditions)
        % mixed block, post-hit and currently non-violation trials
        if contains(conditions{cond}, 'short')
            these = spd.Block==1 & posthit &...
                ~spd.vios & ~postvio & rd_prev<lb;
        elseif contains(conditions{cond}, 'long')
            these = spd.Block==1 & posthit &...
                ~spd.vios & ~postvio & rd_prev>ub;
        end
        rd = nan(size(spd.RewardDelay));
        rd(these) = spd.RewardDelay(these); 
        pokeProb{sess,cond} = nan(length(rd), length(bins));
    
        for ii=1:length(rd)
            poke = [];
            nn = zeros(size(bins));
            if ~isnan(rd(ii))
                % for poke probability
                xx = find(bins<=rd(ii), 1, 'last');
                pokeProb{sess,cond}(ii,1:xx) = 0;
    
                side = spd.RewardedSide{ii}; % get rewarded side
                sOn = speh(ii).States.WaitForSidePoke(1); % get SON time
                if strcmpi(side, 'L')
                    try
                        poke = speh(ii).Events.LeftIn - sOn;
                    catch
                    end
                elseif strcmpi(side, 'R')
                    try
                        poke = speh(ii).Events.RightIn - sOn;
                    catch
                    end
                end
                poke(poke<0) = []; % delete pokes before SON
                poke(poke>=rd(ii))=[]; % delete pokes after SOFF
                % increment N 
                if ~isempty(poke)
                    for jj=1:length(poke)
                        xx = max(find(bins<=poke(jj)));
                        nn(xx) = nn(xx)+1;
                        pokeProb{sess,cond}(ii,xx) = 1;
                    end
                end
                N = N + nn;
    
            end
        end
    end
end

pokeProb_ = {};
for cond=1:length(conditions)
    pokeProb_{cond} = cell2mat(pokeProb(:,cond));
    pokeProb_{cond} = sum(pokeProb_{cond}, 1, 'omitnan') ./ ...
        sum(~isnan(pokeProb_{cond}),1);
end


end



