function ops = glm_defaultOps()
%base function for commonly used forms of all options

ops = struct();

%basics
ops.fname = 'glm_default'; %file save name
ops.dispL = 5; %level of verbosity. 5 gives most output. 1 gives least

%stimulus and data treatment
ops.dt = 0.01; %binsize
ops.T0 = 0; %initial trial time. used only for trial-based structure
ops.T = 1-ops.dt; %final trial time. used only for trial-based structure
ops.usePseudotrial = true; %pseudotrials of time-continuous data, or trial-based structure
ops.stimwrapper = 'gt'; %use the ground-truth form. the most general kind
ops.stimtype = '20200313'; %specific insantiation of ground-truth form
ops.removeData = false; %will remove "opt-out" trials, as defineed by stimtype
ops.nkern = 0; %number of kernels

%basis functions
ops.lh = 0; %length history filter. a 1.2 s filter
ops.ls = ceil(1/ops.dt)*ones(1,ops.nkern); %length stimulus filter. a 1.s s filter.
ops.dim = 0; %dim hist basis
ops.dimst = 7*ones(1,ops.nkern);%dim stim basis. empty. just
ops.dimbc = 0; %dim hist boxcar funcs
ops.dimbcst = 1*ones(1,ops.nkern); %dim stim boxcar funcs
ops.useSpike = 0; %do/do not use spike history kernel
ops.basistype = 'ofc_adaptive'; %use my adaptive basis

%optimization treamtnet
%regularizers. if both 0, will do a grid search in main code to optimize
%them
ops.kfold = 5; %level of cross-validation
ops.cvtype = 'default'; %way to choose cross-val. use trial number
ops.savefreq = 2000; %frequency of temporary saves during opt
ops.funiter = 2000; %max function
ops.gamma = 0*1e-8; %L2 term
ops.beta = 0*1e-6; %L1 term
ops.regularizer = 'L2'; %type of regularizer used
ops.isConstrained = false; %perform a constrained optimization with spike history filter?
ops.useHyperParamGrid = false; %adaptive grid search of the hyperparameters beta and gamma

%stimulus-specific treatments
ops.omitKernels = [];
ops.omitBehavior = [];

%specific options for differnet stim types
ops.spike_wndw = [-ops.T0 0]; %controls binspikes.m for glm data: window (in s) before and after trial to consider spikes. deprecated
ops.ls_supp = 80; %specific to the ofc_adaptive_logUniform
ops.dim_supp = 9; %specific to the ofc_adaptive_logUniform



