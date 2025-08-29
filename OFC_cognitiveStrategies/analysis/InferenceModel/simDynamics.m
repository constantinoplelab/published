function [simControl, simMuscimol] = simDynamics(ratList, behaviorDataPath, ...
    simType, twin, binSize, smoothfactor)
%Simulate wait times and wait time dynamics around block transitions using different
%versions of the behavioral model. SSS 08/2023

%INPUTS:
%   ratList = matfile containing a list of rat names to loop through
%   behaviorDataPath = path to saved behavior data structs
%   simType = string indicating model manipulation to use (e.g. 'lapse')
%   twin = trial window around block transitions
%   binSize = number of trials to bin over. 5 works well
%   smoothfactor = factor for causal smoothing function. 5-7 works well

%OUTPUTS:
%   simControl = simulated control changes in wait-time around block
%       transitions
%   simMuscimol = simulated muscimol changes in wait-time around block
%       transitions

base_params = [0.25 0.30 0.20 0.13]; 
    
for r = 1:length(ratList)
    
    load(strcat(behaviorDataPath, 'ratTrial_', ratList{r},'.mat'));
    A_sim = A;
    A.tau = 2.5;

    if strcmpi(simType, 'lapse')
        lapse = [0.75 0]; %slow updating (muscimol simulation), update on every trial (control simulation)

        for ii = 1:length(lapse)
            [WTOpt, WTMdl] =...
                    GenerateSynthData_Bayes_priorUpdate(base_params, A,...
                    'logn', 1, 8, lapse(ii));  
            A_sim.wait_time = WTOpt;

            [ltom{ii}(r,:), htom{ii}(r,:), mtol{ii}(r,:), mtoh{ii}(r,:), ~, ~, ~] =...
                block_dynamics_wt_binTrials(A_sim, twin, binSize, smoothfactor);
        end

    elseif strcmpi(simType, 'lambda')
        lambdas = [0.25 1]; %sub-optimal prior (muscimol simulation), optimal lambda (control simulation)

        for ii = 1:length(lambdas)
            params = [base_params lambdas(ii)];
            [WTOpt, WTMdl] = ...
                GenerateSynthData_Bayes_subopt_onlylambda_nonflat(...
                params, A, 'logn', 1, 8);

            A_sim.wait_time = WTOpt;

            [ltom{ii}(r,:), htom{ii}(r,:), mtol{ii}(r,:), mtoh{ii}(r,:), ~, ~, ~] =...
                block_dynamics_wt_binTrials(A_sim, twin, binSize, smoothfactor);
        end

    elseif strcmpi(simType, 'kappa')
        kappas = [0.23 0.25 0.22; 0.23 0.30 0.20]; %similar opportunity costs (muscimol simulation), distinct opportunity costs (control simulation)

        for ii = 1:size(kappas, 1)
            params = [kappas(ii,:) 0.13];
            [WTOpt, WTMdl] =...
                GenerateSynthData_Bayes_SS(params, A, 'logn', 1, 8, 'log');

            A_sim.wait_time = WTOpt;

            [ltom{ii}(r,:), htom{ii}(r,:), mtol{ii}(r,:), mtoh{ii}(r,:), ~, ~, ~] =...
                block_dynamics_wt_binTrials(A_sim, twin, binSize, smoothfactor);
        end

    end
end

simControl.mtol = mtol{2};
simControl.ltom = ltom{2};
simControl.mtoh = mtoh{2};
simControl.htom = htom{2};

simMuscimol.mtol = mtol{1};
simMuscimol.ltom = ltom{1};
simMuscimol.mtoh = mtoh{1};
simMuscimol.htom = htom{1};

