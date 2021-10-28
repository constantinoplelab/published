function ops = glm_ofcOps()
%commonly used forms of all options for OFC trials

ops = glm_defaultOps(); %load all of the options

%make the changes for ofc
ops.fname = 'ofc_default';

ops.gamma = 1e-7; %L2 term
ops.beta = 0; %L1 term
%ops.kfoldSupplied = false; %an option to do grid search of hyperparameters. not impelemented in new form yet
ops.removeData = true; %will remove opt-out trials
ops.T = 6; %simulation end time (in s)
ops.T0 = -2; %simulation start time (in s)
ops.dt = 0.05;
ops.usePseudotrial = 1; %use a pure trial-based data strucutre (0, not recommended), or largest contiguous chunks (1, suggested)

%basis function params
ops.cvtype = 'ofc'; %way to choose cross-val. use trial number
ops.basistype = 'ofc_adaptive_logUniform'; %use my adaptive basis
ops.stimtype = '20200710';
ops.stimwrapper = 'ofc';
lh = floor(0.4/ops.dt); % 0.4s filter for spike history
ls =  floor(4/ops.dt); % a 4s filter for stimulu

%these features will be compared and checked for sizing in 
%parse_stimuli_wrapper
ops.nkern = 14; %number of kernels
ops.lh = lh;
ops.ls = ls*ones(1,ops.nkern);
ops.dim = 0; %dim hist basis
ops.dimst = 9*ones(1,ops.nkern); %dim stim basis
ops.dimbc = 0; %dim hist boxcar funcs
ops.dimbcst = 1*ones(1,ops.nkern); %dim stim boxcar funcs
ops.ls_supp = 80; %specific to the ofc_adaptive_logUniform
ops.dim_supp = 9; %specific to the ofc_adaptive_logUniform


ops.dispL = 3; %a verbosity handle for out how output the glm gives




