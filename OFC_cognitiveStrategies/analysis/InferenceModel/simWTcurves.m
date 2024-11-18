function [simC, simM] = simWTcurves(ratTrial, simType)
%Simulate wait times using the behavioral model

%INPUTS:
%   ratTrial = example rat trial data containing reward volumes, trial
%       number, etc
%   simType = string indicating model manipulation to use (e.g. 'lapse')

%OUTPUTS:
%   simC = simulated control wait times
%   simM = simulated muscimol wait times


A_sim = ratTrial;
base_params = [0.23 0.30 0.20 0.13]; 

if strcmpi(simType, 'lapse')
    lapse = [0.75 0]; %slow updating, fast updating
    
    for ii = 1:length(lapse)
        [WTOpt, ~] =...
            GenerateSynthData_Bayes_priorUpdate(base_params, ratTrial,...
            'logn', 1, 8, lapse(ii));
        A_sim.wait_time = WTOpt;

        [hi(ii), lo(ii), ~] = wtcurves(A_sim);

    end

elseif strcmpi(simType, 'lambda')
    lambdas = [0.25 1]; %low lambda, high lambda

    for ii = 1:length(lambdas)
        params = [base_params lambdas(ii)];
        [WTOpt, WTMdl] = ...
            GenerateSynthData_Bayes_subopt_onlylambda_nonflat(...
            params, ratTrial, 'logn', 1, 8);

        A_sim.wait_time = WTOpt;

        [hi(ii), lo(ii), ~] = wtcurves(A_sim);

    end

elseif strcmpi(simType, 'kappa')
    kappas = [0.23 0.25 0.22; 0.23 0.30 0.20];

    for ii = 1:size(kappas,1)
        params = [kappas(ii,:) 0.13];
        [WTOpt, WTMdl] =...
            GenerateSynthData_Bayes(params, ratTrial, 'logn', 1, 8);

        A_sim.wait_time = WTOpt;

        [hi(ii), lo(ii), ~] = wtcurves(A_sim);

    end

end

%wt curves
simC.hi = hi(2).wt;
simC.lo = lo(2).wt;

simM.hi = hi(1).wt;
simM.lo = lo(1).wt;
