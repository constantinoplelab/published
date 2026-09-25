function trials_shuffle = makeshuffle_balanced(Trials,seed)
% Trials
% varargin
    % 1. average across many interations of shuffles
% rng(seed+6)


%combine the two variables
for i = 1:length(Trials)
equate = length(Trials{i,1})-length(Trials{i,2});
if equate>0
removethese = randperm(length(Trials{i,1}));
Trials{i,1}(removethese(1:abs(equate)))=[];
elseif equate<0
removethese = randperm(length(Trials{i,2}));
Trials{i,2}(removethese(1:abs(equate)))=[];  
end
t1 = randperm(length(Trials{i,1}));
t2 = randperm(length(Trials{i,2}));
half = floor(length(t1)/2);
T1 = [Trials{i,1}(t1(1:half)),Trials{i,2}(t2(1:half))];
T1 = reshape(T1,1,numel(T1));
T2 = [Trials{i,1}(t1(half+1:half*2)),Trials{i,2}(t2(half+1:half*2))];
T2 = reshape(T2,1,numel(T2));

Trials_shuffle{i,1} = T1;
Trials_shuffle{i,2} = T2;
end

alignto = {'COFF', 'SON','SOFF','Rew','Opt','CON'};
for k = 1:5
trials_shuffle.(alignto{k}) = Trials_shuffle;
end
