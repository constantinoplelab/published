function [normval, rbar] = divisiveNorm_fun(A, params)
k = params(1);
alpha = params(2);
tback = params(3);


vt = nan(length(A.wait_time), 1);
rbar = vt;
normval = vt;

[con, rew] = convertreward(A.reward);
rew = log2(rew);

for j = 1:length(A.wait_time)
    vt(j) = rew(j);
    if j>tback
        rbar(j) = sum(vt(j-tback:j), 'omitnan');
    else
        rbar(j) = sum(vt(1:j), 'omitnan');
    end
    normval(j) = k*(vt(j)./(1+alpha*rbar(j))); %divisive normalization
end