function plot_fig4h(datadir)
% Plot reward port nose poke probability during the delay period
% during mixed block. N = 16 rats.
%
% INPUTS
% datadir: file path of data downloaded from Zenodo. Link to download can
% be found in README.


load(fullfile(datadir, 'data-published/ratlist.mat'), 'ratList');

pokeProb = cell(length(ratList), 1);

for rr=1:length(ratList)
    rat = ratList{rr};
    fprintf('%i of %i: %s\n', rr, length(ratList), rat)
    
    [pokeProb{rr},bins] = pokePmixed(datadir, rat);
   
end

pokeProb = cell2mat(pokeProb);

figure;
set(gcf, units='inches', position=[9,4,2.2,1.5])
for r=1:length(ratList)
    plot(log2(bins), pokeProb(r,:), color='#a0a0a0', linewidth=0.5)
    hold on
end
xlabel('Time in trial (s, log_{2})')
xlim([1 4])
ylabel('P(rew port poke)')
set(gca, box='off', tickdir='out', fontsize=8, ytick=[0:0.5:1],...
    xtick=[1:1:4])
subtitle(sprintf('N = %i rats', length(ratList)))

% overlay value function
tau = 2.5; 
pcatch = 0.25;
C = 1-pcatch;
R = mean(log2([5,10,20,40,80]));
V = R/tau*C*exp(-bins./tau)./(1-C+C*exp(-bins./tau));
plot(log2(bins), V, 'k', linewidth=1)


end

function [pokeProb, bins] = pokePmixed(datadir, ratname)

file = strcat('SStruct_', ratname, '.mat');
load(fullfile(datadir, 'data-published\S_Structs', file), 'S');
S_final = get_finalTrainingStageSessions(S);

dt = 0.3;
maxT = 100.; 
bins = 0.5:dt:maxT;

N = zeros(size(bins)); % total number of pokes per time bin

for sess=1:length(S_final.pd)
    spd = S_final.pd{sess};
    speh = S_final.peh{sess};

    these = spd.Block==1 & ~spd.vios;
    
    rd = nan(size(spd.RewardDelay));
    rd(these) = spd.RewardDelay(these); 
    pokeProb{sess,1} = nan(length(rd), length(bins));

    for ii=1:length(rd)
        poke = [];
        nn = zeros(size(bins));
        if ~isnan(rd(ii))
            % for poke probability
            xx = find(bins<=rd(ii), 1, 'last');
            pokeProb{sess}(ii,1:xx) = 0;

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
                    pokeProb{sess}(ii,xx) = 1;
                end
            end
            N = N + nn;

        end
    end

end

pokeProb = cell2mat(pokeProb);
pokeProb = sum(pokeProb, 1, 'omitnan') ./ sum(~isnan(pokeProb),1);



end