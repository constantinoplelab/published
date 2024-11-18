function A = splitData(A, nsess)
% Pulls out specific trials from behavior structs
totalT = sum(A.ntrials(1:nsess)); %must start from session 1
bad = find(A.trainingstage(1:totalT) < 9); %remove trials from early training stages before blocks are introduced; only happens in the first session
trials = 1:totalT;
trials(bad) = [];

ns = [1; A.ntrials(1:nsess)];
for ii = 1:nsess   
    numBad = numel(intersect(bad, ns(ii):ns(ii+1)));
    A.ntrials(ii) = A.ntrials(ii) - numBad;
end

A.date = A.date(1:nsess);
A.ntrials = A.ntrials(1:nsess);
A.trainingstage = A.trainingstage(trials);
A.nic = A.nic(trials);
A.catch = A.catch(trials);
A.prob_catch = A.prob_catch(trials);
A.adapt_block = A.adapt_block(trials);
A.test_block = A.test_block(trials);
A.reward = A.reward(trials);
A.reward_delay = A.reward_delay(trials);
A.block = A.block(trials);
A.hits = A.hits(trials);
A.vios = A.vios(trials);
A.optout = A.optout(trials);
A.wait_time = A.wait_time(trials);
A.wait_for_cpoke = A.wait_for_cpoke(trials);
A.zwait_for_cpoke = A.zwait_for_cpoke(trials);
A.trial_num = A.trial_num(trials);
A.side = A.side(trials);
A.lpoke = A.lpoke(trials);
A.rpoke = A.rpoke(trials);
A.cpoke = A.cpoke(trials);
A.lpokedur = A.lpokedur(trials);
A.rpokedur = A.rpokedur(trials);
A.cpokedur = A.cpokedur(trials);
A.rt = A.rt(trials);
A.slrt = A.slrt(trials);
A.ITI = A.ITI(trials);
A.wait_time_unthresholded = A.wait_time_unthresholded(trials);
A.wait_thresh = A.wait_thresh;
A.block_position = A.block_position(trials);