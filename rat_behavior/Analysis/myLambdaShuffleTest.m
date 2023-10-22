function [shuffleDist, trueDiff, pvals, trueParams] =...
    myLambdaShuffleTest(dataLo, dataHi, nPerms, twin)

niters = 10;
xs = -(twin-10):(twin-10);

% myNorm = @(x) x-mean(x(1:twin));
myNorm = @(x) x;

[group1Params, group2Params] = deal(nan(size(dataLo, 1), 4));

dataAll = [dataLo; dataHi];
dataAll = dataAll(:, 11:end-10);
nRats = size(dataAll, 1);

for ii = 1:nPerms
    disp(ii)
    shuffleI = randperm(nRats); % shuffle indicies
    i1 = shuffleI(1:floor(nRats/2)); % first half goes to group 1
    i2 = shuffleI(floor(nRats/2)+1:end); % second half goes to group 2

    group1 = myNorm(mean(dataAll(i1,:)));
    group1Params(ii,:) = my_fit_sigmoid(xs, group1, niters);
    
    group2 = myNorm(mean(dataAll(i2,:)));
    group2Params(ii,:) = my_fit_sigmoid(xs, group2, niters);
end

shuffleDist = group1Params-group2Params;

loParams =...
    my_fit_sigmoid(xs, myNorm(mean(dataLo(:,11:end-10))), niters);
hiParams =...
    my_fit_sigmoid(xs, myNorm(mean(dataHi(:,11:end-10))), niters);

trueParams = [loParams; hiParams];
trueDiff = hiParams-loParams;

pvals = sum(abs(shuffleDist) > abs(trueDiff))/(2*nPerms);
end
