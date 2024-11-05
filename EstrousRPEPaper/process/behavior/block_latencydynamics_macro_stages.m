function Latencydynamics =...
    block_latencydynamics_macro_stages(twin, smoothfactor, prctile_cutoff,...
    f_ratlist, RatBehaviorData)
%% average change in latency around block dynamics. cmc 11/03/21.
%inputs:
%twin = time window, or number of trials around each block transition (40 is good).
%smoothfactor is smoothing factor, because data is noisy. 5 is good.
%minsess = minimum number of sessions in each stage
%wndw = size of moving avg
%threshold = regression slopes need to be above this, was determined from the random shuffle data
%hd = high definition, or the stage that is considered to be high
%estradiol, 1 = proestrus, 2 = proestrus + estrus, 3 = estrus
%outputs: (grouped within cycle groups in Latencydynamics)
%         ltom = mean delta z-scored latency at low to mixed transitions
%         htom = high to mixed
%         mtol = mixed to low
%         mtoh = mixed to high
%For the outputs, each row corresponds to the average delta latency for a
%rat.

%create variables
cycle = {'Proestrus', 'Diestrus'};
frats = length(f_ratlist);
%initialize variables
Latencydynamics = [];
ltom = nan(frats, twin*2+1);
htom = ltom;
mtol = ltom;
mtoh = ltom;

for rat = 1:frats

    ratname = f_ratlist{rat};
    disp([ratname ' ' num2str(rat) ' out of ' num2str(frats)])

    stagedS = RatBehaviorData.(ratname);
    stagedS = stagedS.S;

    for e = 1:length(cycle)

        stageidx = cellfun(@(x) logical(sum(strcmp(x.Stage,...
            cycle{e}))), stagedS.pd); %find where it matches cycle stage(s)
        stageS.pd = stagedS.pd(stageidx);
        stageS.peh = stagedS.peh(stageidx);

        if ~isempty(stageS.pd)

            A = parse_data_from_mysql(stageS, 1, 1);
            A.ITI(A.ITI>prctile(A.ITI, prctile_cutoff)) = NaN;

            [ltom(rat,:), htom(rat,:), mtol(rat,:), mtoh(rat,:)] = ... %decide in here whether or not to z-score
                block_dynamics_latency(A, twin, smoothfactor, 0, 1, 0, 1); %detrend and use causual smoothing

            Latencydynamics.(cycle{e}).ltom(rat,:) = ltom(rat,:);
            Latencydynamics.(cycle{e}).htom(rat,:) = htom(rat,:);
            Latencydynamics.(cycle{e}).mtol(rat,:) = mtol(rat,:);
            Latencydynamics.(cycle{e}).mtoh(rat,:) = mtoh(rat,:);

        else

            Latencydynamics.(cycle{e}).ltom(rat,:) = nan(1,twin*2+1);
            Latencydynamics.(cycle{e}).htom(rat,:) = nan(1,twin*2+1);
            Latencydynamics.(cycle{e}).mtol(rat,:) = nan(1,twin*2+1);
            Latencydynamics.(cycle{e}).mtoh(rat,:) = nan(1,twin*2+1);

        end

    end

end

end