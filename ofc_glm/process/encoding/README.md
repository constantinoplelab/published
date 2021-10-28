# Process functions to run GLM
This is the core code for contructing a GLM model for the ofc dataset. 
 
 ## Contents 

+ addParams.m:     will insert params and update the ops struct
+ add_covariate.m: adds a covariate to an existin stimulus sctruct
+ chooseBasis_ofc.m: a wrapper function around chooseBasis for ofc studies
+ cvalgroup.m: designs labels for cross-validation groups for stratifying indices
+ cvpartition_wrapper.m: decides how to generate labels for cross-validation
+ fun_glmfit.m: core worker function to fit GLMS
+ glm_defaultOps.m: base function for commonly used forms of all options
+ glm_ofcOps.m: commonly used forms of all options for OFC trials
+ loadGlmDat_ofc.m: parses ofc data into Dat structure for glm optimization
+ nll_poiss.m: calcs neg. log-likelihood and gradient, hessian
+ parse_pseudotrials.m: make pseudotrials of continuous-time trial sets
+ parse_stimuli_wrapper.m: decides which parse_stimuli function to use
+ removeOptOut.m: function to remove optouts easily from data, keeping indices of the S
+ removeParams.m:     will remove params and update the ops struct
+ remove_covariate.m: removes a covariate from an existing stimulus struct
+ simulateGLM.m: generate synthetic data from a glm model
+ simulateGLM_logNormal.m: generates synthetic data from a glm, but from lognormal dist
+ stimulusInfo.m:  grabs index and information about a covariate in stim
