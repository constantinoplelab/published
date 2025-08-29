function err = TCA_chooseRank(data, savePath)
% plots the error for a variety of different rank models, used to select
% rank for the model

%% Fit CP Tensor Decomposition

% these commands require that you download Sandia Labs' tensor toolbox:
% https://github.com/andrewssobral/tensor_toolbox/
% https://www.tensortoolbox.org/

% uses ncp from nonnegfac-matlab toolbox by Jinju Kim

% convert data to a tensor object
data = tensor(data);

% find best # of components based on error
n_fits = 10;
err = zeros(n_fits, 20);
test = 1:20;
for r = 1:length(test)
    for ii = 1:n_fits
        [est_factors,~,~] = ncp(data,test(r),'method','hals');

        % store error
        err(ii,r) = norm(full(est_factors) - data)/norm(data);
    end
end

