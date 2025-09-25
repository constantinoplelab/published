function [modelITI, modelValue, ratTrial, modelRPE] = ...
    run_vanillaTD(datadir, modeltype, alpha, scalingfactor)

% Load behavior data of all DA rats
ratList = loadRats(datadir, 'da');

modelITI = struct;
modelValue = struct;
ratTrial = struct;
modelRPE = struct;
for r=1:length(ratList)
    ratname = ratList{r};
    disp(ratname)
    fname = strcat('ratTrial_', ratname, '.mat');
    load(fullfile(datadir, 'data-published', 'A_structs', fname), 'A');
    [modelITI.(ratname), modelValue.(ratname), ratTrial.(ratname),...
        modelRPE.(ratname)] = ...
        vanillaTD_perRat(A, modeltype, alpha, scalingfactor);
end


end