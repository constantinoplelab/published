function [hi, lo, mix, all] = iticurves_SS(A, postUnrewArg)
%Get ITIs for each block. Applies threshold to filter out top 1% of ITIs

l = A.ITI;
postunrew = [0; ~A.hits(1:end-1)]; % post unrewarded trials (Hit = 0)

if postUnrewArg
        % Pull Post-unrewarded Latency
    usethese = ~isnan(l) & postunrew;

else
   usethese = ~isnan(l);
end

L = l(usethese);
L(L>prctile(L,99)) = NaN;
Lz = (L-mean(L, 'omitnan'))./std(L, 'omitnan');

blk = A.block(usethese);


mix.raw = [mean(L(blk==1), 'omitnan'),...
    std(L(blk==1), 'omitnan')./sqrt(sum(blk==3))];
hi.raw = [mean(L(blk==2), 'omitnan'),...
    std(L(blk==2), 'omitnan')./sqrt(sum(blk==2))];
lo.raw = [mean(L(blk==3), 'omitnan'),...
    std(L(blk==3), 'omitnan')./sqrt(sum(blk==3))];
all.raw = [mean(L, 'omitnan'), std(L, 'omitnan')./sqrt(length(L(~isnan(L))))];

mix.z = [mean(Lz(blk==1), 'omitnan'),...
    std(Lz(blk==1), 'omitnan')./sqrt(sum(blk==3))];
hi.z = [mean(Lz(blk==2), 'omitnan'),...
    std(Lz(blk==2), 'omitnan')./sqrt(sum(blk==2))];
lo.z = [mean(Lz(blk==3), 'omitnan'),...
    std(Lz(blk==3), 'omitnan')./sqrt(sum(blk==3))];
all.z= [mean(Lz, 'omitnan'), std(Lz, 'omitnan')./sqrt(length(Lz(~isnan(Lz))))];


end