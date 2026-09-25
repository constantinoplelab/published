function [pval,perlow] = uncertainty_exploratory(all_SU, all_index, all_S)
%% christine's version of MLB
alignto = {'Rew'};
for m = 1:5
    map1{m} = [];
    map2{m} = [];
end
ul = .8; % probability for mix block
ll = .05; %probability for adapt block
params = [1, 1.2, 0.8, .3, 1];
ctr1 = 0;
ctr2 = 0;
for yy = 1:length(all_SU)
    disp(yy)
    sess = all_index{yy,4};
    S = all_S{sess};
    SU = all_SU{yy};
    [rew] = convertreward(S.RewardAmount);

    for m = 1:length(alignto)
        xvec = SU.xvec.(alignto{m});
        hmat = SU.hmat.(alignto{m});
        hmat = (hmat-mean(hmat(:), 'omitnan'))./std(hmat(:), 'omitnan'); %zscore
       
        %add fields that the Bayes script looks for.
        ratTrial = S;
        ratTrial.prob_catch = S.ProbCatch;
        ratTrial.reward = S.RewardAmount;
        ratTrial.ntrials = length(S.ProbCatch);

        [~, ~, ~, Belief, ~, ~] =...
            BayesMdlOpt(params, ratTrial, 'logn', 0, 8);

        pmix = Belief(1,:)';
        phigh = Belief(2,:)';
        plow = Belief(3,:)';
        %mixed block p(mix)> 0.8, p(adapt) > 0.5
        tlow = find(pmix>ul & plow>ll & rew==3); %low
        thigh = find(pmix>ul & phigh>ll & rew==3); %high
        countTlow(yy) = length(tlow)+length(thigh);
        perlow(yy) = countTlow(yy)/ratTrial.ntrials*100;

        for jj = 1:5 % for all rewards
            ix1 = intersect(tlow, find(rew==jj));
            ix2 = intersect(thigh, find(rew==jj));
            % 2nd mostly likely block is low
            if length(ix1)>1 %multiple trials averaged
                map1{jj} = [map1{jj}; mean(hmat(ix1,:), 'omitnan')];
                ctr1 = ctr1+1;
            elseif length(ix1)==1 %single trial saved
                map1{jj} = [map1{jj}; hmat(ix1,:)];
                ctr1 = ctr1+1;
            end
            % 2nd most likely block is high
            if length(ix2)>1
                map2{jj} = [map2{jj}; mean(hmat(ix2,:), 'omitnan')];
                ctr2 = ctr2+1;
            elseif length(ix2)==1
                map2{jj} = [map2{jj}; hmat(ix2,:)];
                ctr2 = ctr2+1;
            end
        end
    end
end

figure;
for mx = 1:2 %high/low
    for jj=3 %mix block
        if mx==1 
            shadedErrorBar(xvec, mean(map1{jj}, 'omitnan'), ...
                std(map1{jj}, 'omitnan')./sqrt(ctr1), 'lineprops', '-b');
        elseif mx==2
            shadedErrorBar(xvec, mean(map2{jj}, 'omitnan'), ...
                std(map2{jj}, 'omitnan')./sqrt(ctr2), 'lineprops', '-r');
        end
        title(strcat(['p(mix)>', num2str(ul), 'p(adapt)>', num2str(ll)]));
    end

end
    
    for i = 1:length(xvec)
        [pval(i)] = permute_test(map1{3}(:,i),map2{3}(:,i),1000);
    end

    hold on;
    sigT = find(pval<0.05);
    sigd = ones(1,length(xvec));
    plot(xvec(sigT),sigd(sigT)*0.5,'*k');
    % xline(xvec(sigT))
    
    xline(0, '--k')

    %edit figure to be pretty
set(gca, 'TickDir', 'out'); box off
    xlabel('Time from reward (s)');
    ylabel('Firing rate (Hz)');
    xlim([-1 2.5]);xticks(-1:2)
    ylim([-.25 .55]);
    lgd = legend('low block','high block');
    title(lgd, '2nd most likely block');
    set(gcf, 'Color', [1 1 1]);

    
