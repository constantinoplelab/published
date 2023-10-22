function [mean_iti, sem_iti, p] = iti_by_reward(data)

    first_trial = [1; cumsum(data.ntrials(1:end-1))+1];
    data.ITI(first_trial) = nan;
    data.ITI(data.ITI>prctile(data.ITI, 99)) = nan;
    
    post_vios = logical([0; data.vios(1:end-1)]);
    usethese = ~post_vios & data.block==1;

    l = data.ITI(usethese);
    l = (l-mean(l, 'omitnan'))./std(l, 'omitnan');
    
    [~, rew] = convertreward(data.reward(usethese));
    blk = data.block(usethese);
    
    rvec = [5 10 20 40 80];
    mean_iti = arrayfun(@(r) mean(l(rew==r & blk==1), 'omitnan'), rvec);
    sem_iti = arrayfun(@(r)...
        std(l(rew==r)./sqrt(sum(rew==r & blk==1)), 'omitnan'), rvec);
    p = anova1(l(blk==1), rew(blk==1), 'off');
end