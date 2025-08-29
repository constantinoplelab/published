function [d, d_shuff] = compute_dPrime(SU, alignto, wndw, trialtype1, trialtype2)
% Compute the dscriminability index (d') for trialtype1 versus trialtype 2.
% Shuffled d' is subtracted

% INPUTS: 
%   SU = struct of spiking info for a single recording session
%   alignto = task event ('COFF' = trial start, 'SON' = LED on, 'SOFF' =
%       LED off, 'Rew' = reward)
%   wndw = time window around event to use for d' computation
%   trialtype1 = vector of trial indices for one condition (e.g. small
%       volume trials)
%   trialtype2 = vector of trial indices for the second condition (e.g.
%       large volume trials)

% d' formula
dPrime = @(mu1,mu2,var1,var2) abs((mu1 - mu2))./...
    sqrt(.5 * (var1 + var2));

w = arrayfun(@(x) find(round(SU(1).xvec.COFF, 3) == x), wndw);
n = length(SU);

% balance trials
min_t = min([length(trialtype1), length(trialtype2)]);
trialtype1 = trialtype1(randsample(length(trialtype1),min_t));
trialtype2 = trialtype2(randsample(length(trialtype2),min_t));

% pull out binned firing rates aligned to task event of interest
hmat1 = arrayfun(@(x) SU(x).hmat.(alignto)(trialtype1, w(1):w(2)), 1:n, ...
    'UniformOutput', false);
hmat2 = arrayfun(@(x) SU(x).hmat.(alignto)(trialtype2, w(1):w(2)), 1:n, ...
    'UniformOutput', false);

mu1 = cell2mat(arrayfun(@(x) mean(hmat1{x}, 'omitnan'), 1:n, 'UniformOutput', false)');
var1 = cell2mat(arrayfun(@(x) var(hmat1{x}, 0, 'omitnan'), 1:n, 'uniformoutput', false)');

mu2 = cell2mat(arrayfun(@(x) mean(hmat2{x}, 'omitnan'), 1:n, 'UniformOutput', false)');
var2 = cell2mat(arrayfun(@(x) var(hmat2{x}, 0, 'omitnan'), 1:n, 'uniformoutput', false)');

dprime = dPrime(mu1, mu2, var1, var2);

% shuffles
HMAT = arrayfun(@(x) [hmat1{x}; hmat2{x}], 1:n, 'UniformOutput', false);
a = size(HMAT{1}, 1); %all should be the same size

nshuff = 15;
dshf = nan(n, length(w(1):w(2)), nshuff); %cells x time x shuffles

for j = 1:nshuff
    these1 = randperm(a, min_t);
    these2 = randperm(a, min_t);

    mu1_shuff = cell2mat(arrayfun(@(x) mean(HMAT{x}(these1,:),'omitnan'), 1:n, ...
        'uniformoutput', false)');
    var1_shuff = cell2mat(arrayfun(@(x) var(HMAT{x}(these1,:), 0,'omitnan'), 1:n, ...
        'uniformoutput', false)');

    mu2_shuff = cell2mat(arrayfun(@(x) mean(HMAT{x}(these2,:),'omitnan'), 1:n, ...
        'uniformoutput', false)');
    var2_shuff = cell2mat(arrayfun(@(x) var(HMAT{x}(these2,:), 0,'omitnan'), 1:n, ...
        'UniformOutput', false)');

    dshf(:,:,j) = dPrime(mu1_shuff, mu2_shuff, var1_shuff, var2_shuff);
end

%subtract shuffled d'   
d_shuff = mean(dshf, 3, 'omitnan');
d = dprime - d_shuff;

end

