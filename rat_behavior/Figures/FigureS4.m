function FigureS4(datadir, codedir)
%FigureS4 - Average trial initiation time in mixed blocks conditioned on
% the previous block
%   datadir = directory of dataset
%   codedir = director of code

%% Set up paths
s = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr,...
    codedir, 'IgnoreCase', ispc);

if ~onPath % only add code dir to path if it already isn't
    addpath(genpath(codedir))
end

%% Process data

a = load([datadir filesep 'ratList.mat']);
ratList = a.ratList; % list of rats

% Pre-allocate matricies
MixedITIbyPrevBlock = nan(length(ratList), 2);

for rr = 1:length(ratList)
    fprintf('%d out of %d\n', rr, length(ratList));

    % Load data
    A = load([datadir 'A_Structs_Final' filesep...
        'ratTrial_' ratList{rr} '.mat']);
    A = A.A;

    % Pull and z-score latency
    L = A.ITI;
    L(L>prctile(L, 99)) = nan;
    L = (L-mean(L, 'omitnan'))./std(L, 'omitnan');

    % Find previous block
    blockTrans = [nan; diff(A.block)]; % find all block transitions
    mixedBlockTrans = find(blockTrans<0); % block transitions into mixed
    mixedBlockI = find(A.block==1); % convert to indicies

    [PrevLoBloc, PrevHiBloc] = deal({}); % Preallocate

    % Loop over each block transition
    for ii = 1:length(mixedBlockTrans)-1
        usethese = intersect(mixedBlockTrans(ii):...
            mixedBlockTrans(ii+1)-1,...
            mixedBlockI);

        if blockTrans(mixedBlockTrans(ii)) == -1
            PrevHiBloc{end+1} = L(usethese); %#ok<*SAGROW>
        elseif  blockTrans(mixedBlockTrans(ii)) == -2
            PrevLoBloc{end+1} = L(usethese);
        end
    end
    PrevHiBloc = cell2mat(PrevHiBloc'); % Find mixed blks after high blks
    PrevLoBloc = cell2mat(PrevLoBloc'); % Find mixed blks after low blks

    % Average latency
    MixedITIbyPrevBlock(rr,:) =...
        [mean(PrevHiBloc, 'omitnan'), mean(PrevLoBloc, 'omitnan')];
end

%% Plot 

figure; hold on
plot(1:2, MixedITIbyPrevBlock, color=[0.7 0.7 0.7])
plot(1:2, mean(MixedITIbyPrevBlock), 'k')

xlim([0.5 2.5])
ylim([-0.23 0.23])

xticks(1:2)
xticklabels({'High', 'Low'})

xlabel('Previous Block')
ylabel('Trial initiation time (z-score)')

title('Mixed Block ITI')

shg
end