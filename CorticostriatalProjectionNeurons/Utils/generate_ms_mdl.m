function [itiMdl, itiOpt, V, rpe, alpha, G] = generate_ms_mdl(...
    R, alpha0, D)
[V, rpe, alpha] = deal(nan(size(R)));
V(1) = 2.5;

alphas = alpha0*(1:5);

for tt = 1:length(R)
    alpha(tt) = alphas(R(tt));

    rpe(tt) = R(tt) - V(tt);
    V(tt+1) = V(tt) + alpha(tt)*rpe(tt);
end
G = alpha./alpha0;

itiOpt = D./V;
itiMdl = generate_wts(itiOpt, 1, 'logn', 1.5);

end