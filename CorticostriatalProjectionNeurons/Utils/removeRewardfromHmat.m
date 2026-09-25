function all_SU_removed = removeRewardfromHmat(goods, all_SU,all_S,all_index)

xvec = all_SU{goods(1)}.xvec.COFF;
all_SU_removed = all_SU;
%nan things that are after the next event
for a= 1:length(goods)
S = all_S{all_index{goods(a),4}};
SU = all_SU{goods(a)};
%find delay
wt = S.wait_time;
wt(isnan(wt)) = 0;
wtT = [];
wtTCOff = [];
for t = 1:length(wt)
wtT(t,:) = xvec>=wt(t);
wtTCOff(t,:) = xvec>=(wt(t)+S.NoseInCenter(t));
end
%Nan anything after the delay
hmat = SU.hmat.SON;
hmat(find(wtT)) = NaN;
all_SU_removed{goods(a)}.hmat.SON = hmat;

hmat = SU.hmat.COFF;
hmat(find(wtTCOff)) = NaN;
all_SU_removed{goods(a)}.hmat.COFF = hmat;
end
