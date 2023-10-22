function [hi, lo, mix, p] = iticurves(A)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

l = A.ITI;
usethese = ~isnan(l);

L = l(usethese);
L(L>prctile(L,99)) = NaN;
Lz = (L-mean(L, 'omitnan'))./std(L, 'omitnan');

blk = A.block(usethese);
p = ranksum(L(blk==2), L(blk==3), tail='left');

mix.raw = [mean(L(blk==1), 'omitnan'),...
    std(L(blk==1), 'omitnan')./sqrt(sum(blk==3))];
hi.raw = [mean(L(blk==2), 'omitnan'),...
    std(L(blk==2), 'omitnan')./sqrt(sum(blk==2))];
lo.raw = [mean(L(blk==3), 'omitnan'),...
    std(L(blk==3), 'omitnan')./sqrt(sum(blk==3))];

mix.z = [mean(Lz(blk==1), 'omitnan'),...
    std(Lz(blk==1), 'omitnan')./sqrt(sum(blk==3))];
hi.z = [mean(Lz(blk==2), 'omitnan'),...
    std(Lz(blk==2), 'omitnan')./sqrt(sum(blk==2))];
lo.z = [mean(Lz(blk==3), 'omitnan'),...
    std(Lz(blk==3), 'omitnan')./sqrt(sum(blk==3))];

end