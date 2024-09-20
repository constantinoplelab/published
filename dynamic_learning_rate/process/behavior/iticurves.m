function [hi, lo, mix, p] = iticurves(A)
%iticurves - Average trial intiation time by block
% INPUTS:
%   A - Rat behavioral data. Usually BestFitEarly.(RAT_ID).All.ratTrial
% OUTPUTS:
%   hi, lo, mix - First row is average ITI in second (raw) and z-scored (z) 
%       for each rat in high block (hi), low block (lo), or mixed block 
%       (mix). Second row is SEM.
%   p - p-value for pairwise comparisons. Order is High vs. Low, Mixed vs.
%   High, and Mixed vs. Low

% Pull ITIs
l = A.ITI;

usethese = ~isnan(l); % Remove nan
L = l(usethese);

L(L>prctile(L,99)) = NaN; % Remove outliers (99th percentile)

Lz = (L - mean(L, 'omitnan')) ./ std(L, 'omitnan'); % Z-score

blk = A.block(usethese); % Block vector

% Do pairwise rank-sum tests
p = [ranksum(L(blk==2), L(blk==3)),...
    ranksum(L(blk==1), L(blk==2)),...
    ranksum(L(blk==1), L(blk==3))];

% Average raw ITI
mix.raw = [mean(L(blk==1), 'omitnan'),...
    std(L(blk==1), 'omitnan')./sqrt(sum(blk==1))];

hi.raw = [mean(L(blk==2), 'omitnan'),...
    std(L(blk==2), 'omitnan')./sqrt(sum(blk==2))];

lo.raw = [mean(L(blk==3), 'omitnan'),...
    std(L(blk==3), 'omitnan')./sqrt(sum(blk==3))];

% Average z-score ITI
mix.z = [mean(Lz(blk==1), 'omitnan'),...
    std(Lz(blk==1), 'omitnan')./sqrt(sum(blk==1))];

hi.z = [mean(Lz(blk==2), 'omitnan'),...
    std(Lz(blk==2), 'omitnan')./sqrt(sum(blk==2))];

lo.z = [mean(Lz(blk==3), 'omitnan'),...
    std(Lz(blk==3), 'omitnan')./sqrt(sum(blk==3))];
