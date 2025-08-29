function A = splitData(A, nsess)
% Pulls out specific trials from behavior structs
sessions = 1:nsess;
totalT = sum(A.ntrials(1:nsess));
bad = 1:find(A.trainingstage(1:totalT) < 9, 1, 'last'); %remove trials from early training stages before blocks are introduced; only happens in the first session

ns = cumsum([1; A.ntrials(sessions)]);
for ii = 1:nsess   
    numBad = numel(intersect(bad, [ns(ii):ns(ii+1)-1]));
    ntrials(ii) = A.ntrials(ii) - numBad;
end

iszero = find(ntrials(sessions) == 0);
sessions(iszero) = []; %get rid of any sessions without any stage 9 sessions (only happens if the rat is brought back down a stage)
if ~isempty(iszero)
    sessions = sessions(1):sessions(1) + nsess-1;
end

totalT = sum(A.ntrials(1:sessions(end))); %must start from session 1, need to index from original trials
trials = 1:totalT;
trials(bad) = [];
A.ntrials(1:nsess) = ntrials;

A.date = A.date(sessions);
A.ntrials = A.ntrials(sessions);
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
A.trial_num = A.trial_num(trials);
A.side = A.side(trials);
A.rt = A.rt(trials);
A.ITI = A.ITI(trials);
try
    A.block_position = A.block_position(trials);
    A.zwait_for_cpoke = A.zwait_for_cpoke(trials);
    A.lpoke = A.lpoke(trials);
    A.rpoke = A.rpoke(trials);
    A.cpoke = A.cpoke(trials);
    A.lpokedur = A.lpokedur(trials);
    A.rpokedur = A.rpokedur(trials);
    A.cpokedur = A.cpokedur(trials);
    A.slrt = A.slrt(trials);
    A.wait_time_unthresholded = A.wait_time_unthresholded(trials);
catch
end
A.wait_thresh = A.wait_thresh;
    
end