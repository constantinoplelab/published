function [tensor_idx, mdl_idx] = findZeroNF(model, IDmap)
%finds cells with a neuron factor of 0

%INPUTS:
%   model = ktensor, this is already sorted within components
%   IDmap = provides a mapping of the neurons in the model

%OUTPUTS:
%   tensor_idx = index of the cells with 0 nf in the tensor as constructed
%       and cellInfo

nf = model.U{1}; %neuron factors

mdl_idx = find(sum(nf,2) < 0.0001 );

tensor_idx = IDmap(mdl_idx,1);





