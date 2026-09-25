function [itiMdl, itiOpt, V, rpe, alpha, G] = generate_ph_mdl(...
    R, alpha0, D)
[V, rpe, alpha, G] = deal(nan(size(R)));
V(1) = 2.5;

for tt = 1:length(R)
    if tt == 1
        G(tt) = 1;
    else
        G(tt) = abs(rpe(tt-1));
    end
    alpha(tt) = min(alpha0*G(tt), 1);

    rpe(tt) = R(tt) - V(tt);
    V(tt+1) = V(tt) + alpha(tt)*rpe(tt);
end

itiOpt = D./V;
itiMdl = generate_wts(itiOpt, 1, 'logn', 1.5);

end