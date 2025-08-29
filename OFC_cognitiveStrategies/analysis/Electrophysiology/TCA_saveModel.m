function TCA_saveModel(tensorPath, savePath, R, titles)
%save the model with the lowest error and corresponding figure
%INPUTS:
%   tensorPath = path to tensors 
%   savePath = path to save in
%   R = number of components
%   titles = plotting titles e.g. {'neurons' 'time' 'trials'};

matfiles = dir(fullfile(tensorPath,'*Tensor.mat'));
info = dir(fullfile(tensorPath, '*Info.mat'));

% xvec = repmat(-0.5:.05:1, 1, 5);
x2 = repmat(-1:.05:0.8, 1, 2);
xvec = [repmat(-1:.05:1, 1, 3) x2];
drawLines = nan(2,5);
drawLines(1,:) = find(xvec==0);
drawLines(2,1:3) = find(xvec == 1);
drawLines(2,4:5) = find(xvec == 0.8, 2, 'last');

for s = 1:length(matfiles)
    
    data = load([matfiles(s).folder filesep matfiles(s).name]);
    info = load([info(s).folder filesep info(s).name]);

    name = extractBefore(matfiles(s).name,'.');
    f = fieldnames(data);
    fi = fieldnames(info);
    % convert data to a tensor object
    data = tensor(data.(f{1}));
    info = info.(fi{1});

    [model,IDmap] = runTCA(data,R,titles,0);

    %add tca cluster ID to info
    ncomp = max(IDmap(:,2));

    info.TCAclust = nan(length(IDmap),1);
    info.TCAclust(IDmap(:,1)) = IDmap(:,2);

    %relabel neurons with neuron factors of 0 for all components
    [tensor_idx, mdl_idx] = findZeroNF(model, IDmap);
    info.TCAclust(tensor_idx) = 0;

    % if size(data,3) > 3
    %     plotTCA(model, drawLines,  ...
    %         'Plottype', {'bar', 'line', 'line'}, ...
    %         'Modetitles', titles);
    % else
    %     plotTCA(model, drawLines,  ...
    %         'Plottype', {'bar', 'line', 'line'}, ...
    %         'Modetitles', titles);
    % end
    % set(gcf, 'Name', ['estimated factors - fit #' num2str(n)])
    % set(gcf,'renderer','painter')
    
    % save(strcat(savePath,'model_',num2str(R),'comp_',matfiles(s).name),...
    %     'model')
    % save(strcat(savePath,'IDmap_',num2str(R),'comp'),'IDmap')
    % save(strcat(savePath, 'info_',num2str(R),'comp'), 'info')
    % saveas(gcf,strcat(savePath,num2str(R),'comp_',name,'.fig'))
end