function [basefun,basefunst] = chooseBasis_ofc(ops)
%CHOOSEBASIS_OFC: a wrapper function around chooseBasis for ofc studies

%length of history filter in timesteps, 
lh = ops.lh; %length history filter. a 1.2 s filter
ls = ops.ls; %length stimulus filter. a 1.s s filter.
dim = ops.dim; %dim hist basis
dimst = ops.dimst; %dim stim basis
dimbc = ops.dimbc; %dim hist boxcar funcs
dimbcst = ops.dimbcst; %dim stim boxcar funcs


if ~isfield(ops,'basistype')
    basistype = 'ofc_adaptive';
else
    basistype = ops.basistype;
end

%sopike history basis basis function
%TODO: REMOVE ofc_adaptive plus stuff
if strcmp(basistype,'ofc_adaptive_logUniform') || strcmp(basistype,'ofc_adaptive_plus')
    bf_params_hist = [lh,dim,dimbc,ops.ls_supp,ops.dim_supp];
else
    bf_params_hist = [lh,dim,dimbc];
end
basefun = chooseBasis(basistype,1:ops.lh,bf_params_hist);

%some contingencies for if spike history kernel was not wanted, but basis
%functions were accidentally defined with number basis funs, etc.
if ~ops.useSpike
    %basefun = 0*basefun;
    basefun = [];
else
    basefun(:,1) = 0; %make sure reverse time form has 0 at initial pt.
end


%stimulus basis functions
%TODO: change to number of kernels, not covariates

%basefunst = cell(ops.ncov-1,1);
basefunst = cell(ops.nkern,1);
for j = 1:ops.nkern
    %TODO: REMOVE ofc_adaptive plus stuff
    if strcmp(basistype,'ofc_adaptive_logUniform')  || strcmp(basistype,'ofc_adaptive_plus')%requires a reference ls to place basis functions
        bf_params_stim = [ls(j),dimst(j),dimbcst(j),ops.ls_supp,ops.dim_supp];
        basefunst{j} = chooseBasis(basistype,1:ops.ls_supp,bf_params_stim);
    elseif numel(ls)==1 && numel(dimst)==1 && numel(dimbcst)==1 %only one set of info provided for all basis functions. all same
        bf_params_stim = [ls,dimst,dimbcst];
        basefunst{j} = chooseBasis(basistype,1:ls(j),bf_params_stim);
    else
        bf_params_stim = [ls(j),dimst(j),dimbcst(j)];
        basefunst{j} = chooseBasis(basistype,1:ls(j),bf_params_stim);
    end
    
    %basefunst{j} = chooseBasis(basistype,1:ls(j),bf_params_stim);
end


