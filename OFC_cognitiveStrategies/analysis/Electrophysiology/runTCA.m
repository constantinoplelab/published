function [model,IDmap] = runTCA(data,R,titles,plotArg)
% converts data matrix to tensor object then finds the model with the
% lowest error for the provided rank. Plots model

%Does not remove neurons with neuron factors of 0 at this stage

%INPUTS:
%   data = data matrix created in time__Tensor
%   R = number of components to use in model determined by TCA_chooseRank
%   titles = plotting titles e.g. {'neurons' 'time' 'blocks'}
%   plotArg = 1 to plot

%OUPUTS:
%   model: ktensor
%   IDmap: mapping from original data matrix/order in cellInfo table to
%       model order

%modified by SSS 2022
%% Fit CP Tensor Decomposition

% these commands require that you download Sandia Labs' tensor toolbox:
% https://github.com/andrewssobral/tensor_toolbox/
% https://www.tensortoolbox.org/

com = @(m,x) sum(m.*x) / sum(m);
% convert data to a tensor object
data = tensor(data);
     
%similarity
% score(A,B,varargin)
xvec = repmat(-0.5:.05:1, 1, 5);
drawLines(1,:) = find(xvec==0);
drawLines(2,:) = find(xvec == 1);

n_fits = 10;
% for R = 1:7
for n = 1:n_fits
    % fit model
    [model{n},~,info{n}] = ncp(data,R,'method','hals'); % non-negative decomp
end

%select model with the lowest reconstruction error
modelEr = arrayfun(@(x) info{x}.final.rel_Error, 1:length(info));
minEr_idx = find(modelEr == min(modelEr), 1);
model = model{minEr_idx};

%sort time factors by time to center of mass, rearrange all orders to match
%neurons stay in the same order but their weights are ordered by center of
%mass
comByComponent = arrayfun(@(x) com(model.U{2}(:,x),[1:length(model.U{2})]'), ...
    1:R);
[~,order] = sort(comByComponent,'ascend');
for ii = 1:3
    model.U{ii} = model.U{ii}(:,order);
end

%then sort neuron factors by their max component weight
[sorted, IDmap] = sortNeuronFact(model,R); %IDmap keeps track of the original order that is in cellInfo 
model.U{1} = sorted;

if plotArg
    plotTCA(model, drawLines,  ...
        'Plottype', {'bar', 'line', 'line'}, ...
        'Modetitles', titles);
    set(gcf, 'Name', ['estimated factors - fit #' num2str(n)])
    set(gcf,'renderer','painter')
end


disp(min(modelEr))