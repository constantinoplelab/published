function [data, beta, popt] = parse_data_from_mysql(S, varargin)
%concatenate data over days.
data.date = [];
data.trainingstage = [];
data.nic = [];
data.catch = [];
data.prob_catch = [];
data.adapt_block = [];
data.test_block = [];
data.reward = [];
data.reward_delay = [];
data.block = [];
data.hits = [];
data.vios = [];
data.optout = [];
data.wait_time = [];
data.wait_for_cpoke = [];
data.zwait_for_cpoke = []; %z-scored latency
data.timeout = [];
data.ntrials = [];
data.trial_num = [];
data.side = [];
data.lpoke = [];
data.rpoke = [];
data.cpoke = [];
data.lpokedur = [];
data.rpokedur = [];
data.cpokedur = [];
data.rt = [];
data.slrt = []; % side led reaction time; added by HJ 2022/11
data.firstpoke = [];

struct_cond = cellfun(@isstruct, S.pd);
S.pd = S.pd(struct_cond);
S.peh = S.peh(struct_cond);

% Check screen perf unless you just want to look at one date
if isempty(varargin)
    [~, beta, goods] = screen_perf(S);
    [popt] = screen_optrate(S);
elseif length(varargin)==2
    goods = 1:length(S.pd);
end
    
if ~isempty(varargin)
    if ~isempty(varargin{1})
        beta = ones(length(S.pd), 1);
    end
end

prot_cond = cellfun(@(s) ~strcmp(s.Protocol, 'RewardWaitTime'), S.pd);
hit_cond = cellfun(@(s) length(s.hits)>=100, S.pd);
has_data = cellfun(@(s) isfield(s, 'wait_time'), S.pd);
all_conds = find(prot_cond & hit_cond & has_data);
include_sess = intersect(all_conds, goods);

for ii = 1:length(include_sess)
    j = include_sess(ii);
    data.date = [data.date;...
        datetime(S.pd{j}.SessionDate, 'InputFormat', 'dd-MMM-yyyy')];
    data.trainingstage = [data.trainingstage; S.pd{j}.TrainingStage];
    data.nic = [data.nic; S.pd{j}.NoseInCenter];
    data.catch = [data.catch; S.pd{j}.RewardDelay==100];
    data.prob_catch = [data.prob_catch; S.pd{j}.ProbCatch];
    data.adapt_block = [data.adapt_block; S.pd{j}.BlockLengthAd];
    data.test_block = [data.test_block; S.pd{j}.BlockLengthTest];
    data.reward = [data.reward; S.pd{j}.RewardAmount];
    data.reward_delay = [data.reward_delay; S.pd{j}.RewardDelay];
    data.block = [data.block; S.pd{j}.Block];
    data.hits = [data.hits; S.pd{j}.hits];
    data.vios = [data.vios; S.pd{j}.vios];
    data.optout = [data.optout; S.pd{j}.optout];
    data.wait_time = [data.wait_time; S.pd{j}.wait_time];
    data.wait_for_cpoke = [data.wait_for_cpoke;...
        nan; S.pd{j}.WaitForPoke(2:end)];
    data.zwait_for_cpoke = [data.zwait_for_cpoke; nan; ...
        (S.pd{j}.WaitForPoke(2:end) -...
        nanmean(S.pd{j}.WaitForPoke(2:end)))./...
        nanstd(S.pd{j}.WaitForPoke(2:end))];
    data.ntrials = [data.ntrials; length(S.pd{j}.hits)];
    data.trial_num = [data.trial_num; (1:length(S.pd{j}.hits))'];
    data.rt = [data.rt; S.pd{j}.ReactionTime];
    
    for mm = 1:length(S.pd{j}.hits)
        if  ~isempty(S.peh)
            if isfield(S.peh{j}(mm).Events, 'CenterIn') &&...
                    isfield(S.peh{j}(mm).Events, 'CenterOut')
                if length(S.peh{j}(mm).Events.CenterOut)<...
                        length(S.peh{j}(mm).Events.CenterIn)
                    if S.peh{j}(mm).Events.CenterIn(end)>...
                            S.peh{j}(mm).Events.CenterOut(end)
                        S.peh{j}(mm).Events.CenterIn(end) = [];
                    end
                elseif length(S.peh{j}(mm).Events.CenterOut)>...
                        length(S.peh{j}(mm).Events.CenterIn)
                    if S.peh{j}(mm).Events.CenterIn(1)>...
                            S.peh{j}(mm).Events.CenterOut(1)
                        S.peh{j}(mm).Events.CenterOut(1) = [];
                    end
                end
                data.cpokedur = [data.cpokedur; ...
                    sum(S.peh{j}(mm).Events.CenterOut(2:end)-...
                    S.peh{j}(mm).Events.CenterIn(2:end))];
                data.cpoke = [data.cpoke;...
                    length(S.peh{j}(mm).Events.CenterOut)];
                
                if isfield(S.peh{j}(mm).Events, 'LeftIn') &&...
                        isfield(S.peh{j}(mm).Events, 'LeftOut')
                    if length(S.peh{j}(mm).Events.LeftOut)<...
                            length(S.peh{j}(mm).Events.LeftIn)
                        if S.peh{j}(mm).Events.LeftIn(end)>...
                                S.peh{j}(mm).Events.LeftOut(end)
                            S.peh{j}(mm).Events.LeftIn(end) = [];
                        end
                    elseif length(S.peh{j}(mm).Events.LeftOut)>...
                            length(S.peh{j}(mm).Events.LeftIn)
                        if S.peh{j}(mm).Events.LeftIn(1)>...
                                S.peh{j}(mm).Events.LeftOut(1)
                            S.peh{j}(mm).Events.LeftOut(1) = [];
                        end
                    end
                    usethese = find(S.peh{j}(mm).Events.LeftIn>...
                        S.peh{j}(mm).Events.CenterOut(1));
                    data.lpokedur = [data.lpokedur; ...
                        sum(S.peh{j}(mm).Events.LeftOut(usethese)-...
                        S.peh{j}(mm).Events.LeftIn(usethese))];
                    data.lpoke = [data.lpoke;...
                        length(S.peh{j}(mm).Events.LeftOut(usethese))];
                else
                    data.lpokedur = [data.lpokedur; nan];
                    data.lpoke = [data.lpoke; nan];
                end
                
                %  if isfield(S.peh{j}(mm).Events, 'CenterIn') &&...
                %isfield(S.peh{j}(mm).Events, 'CenterOut')
                
                if isfield(S.peh{j}(mm).Events, 'RightIn') &&...
                        isfield(S.peh{j}(mm).Events, 'RightOut')
                    if length(S.peh{j}(mm).Events.RightOut) <...
                            length(S.peh{j}(mm).Events.RightIn)
                        if S.peh{j}(mm).Events.RightIn(end)>...
                                S.peh{j}(mm).Events.RightOut(end)
                            S.peh{j}(mm).Events.RightIn(end) = [];
                        end
                    elseif length(S.peh{j}(mm).Events.RightOut) >...
                            length(S.peh{j}(mm).Events.RightIn);
                        if S.peh{j}(mm).Events.RightIn(1)>...
                                S.peh{j}(mm).Events.RightOut(1)
                            S.peh{j}(mm).Events.RightOut(1) = [];
                        end
                    end
                    usethese = find(S.peh{j}(mm).Events.RightIn>...
                        S.peh{j}(mm).Events.CenterOut(1));
                    data.rpokedur = [data.rpokedur; ...
                        sum(S.peh{j}(mm).Events.RightOut(usethese) -...
                        S.peh{j}(mm).Events.RightIn(usethese))];
                    data.rpoke = [data.rpoke;...
                        length(S.peh{j}(mm).Events.RightOut(usethese))];
                    
                else
                    data.rpokedur = [data.rpokedur; nan];
                    data.rpoke = [data.rpoke; nan];
                    
                end
            else
                data.cpokedur = [data.cpokedur; nan];
                data.cpoke = [data.cpoke; nan];
                data.lpokedur = [data.lpokedur; nan];
                data.lpoke = [data.lpoke; nan];
                data.rpokedur = [data.rpokedur; nan];
                data.rpoke = [data.rpoke; nan];
            end
        else
            data.cpokedur = [data.cpokedur; nan];
            data.cpoke = [data.cpoke; nan];
            data.lpokedur = [data.lpokedur; nan];
            data.lpoke = [data.lpoke; nan];
            data.rpokedur = [data.rpokedur; nan];
            data.rpoke = [data.rpoke; nan];
        end
        data.side = [data.side; S.pd{j}.RewardedSide{mm}];
    end
end

% problem_fields = {'trainingstage', 'nic', 'catch', 'prob_catch',...
%     'adapt_block', 'test_block', 'reward', 'reward_delay', 'block'};
% for ff = 1:length(problem_fields)
%    data.(problem_fields{ff})(isnan(data.(problem_fields{ff}))) = [];
% end

ITI = getITI(structfun(@(s) s(include_sess), S, 'UniformOutput', false));
data.ITI = cell2mat(ITI');

slrt = getSLRT(structfun(@(s) s(include_sess), S, 'UniformOutput', false));
data.slrt = cell2mat(slrt');

%if the trial timed out, we don't want to treat it like a true catch trial.
data.catch(data.timeout==1) = nan;
data.reward_delay(data.timeout==1) = nan;

data.wait_time_unthresholded = data.wait_time;

%exclude outlier wait times greater than 1 SD above the mean.
data.wait_thresh = mean(data.wait_time(data.reward_delay==100),...
    'omitnan') + std(data.wait_time(data.reward_delay==100), 'omitnan');
data.wait_time(data.wait_time>data.wait_thresh) = nan;
end
