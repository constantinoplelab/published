function [itiMdl, itiOpt, V, rpe, alpha, G, belief] =...
    generate_deltabelief_mdl(...
    R, alpha0, D)

belief = nan(3, length(R));
belief(:,1) = [1 0 0];

p_rew = nan(3, 5);
p_rew(1,:) = [1/5 1/5 1/5 1/5 1/5]; % mixed block
p_rew(2,:) = [0 0 1/3 1/3 1/3]; % high block
p_rew(3,:) = [1/3 1/3 1/3 0 0]; % low block

H0 = 1/40;
p_b = [1-H0, H0, H0;...
    H0/2, 1-H0, 0;...
    H0/2 0, 1-H0];

for rr = 2:length(R)
    likey = p_rew(:, R(rr));
    prior = p_b*belief(:, rr-1);

    belief(:,rr) = likey.*prior;
    belief(:,rr) = belief(:,rr)./sum(belief(:,rr));
end

G = 1./(1-abs(belief(1,1:end-1) - belief(1,2:end)));

[V, alpha, rpe] = deal(nan(size(R)));
V(1) = 2.5;

for tt = 1:length(R)-1
    rpe(tt) = R(tt) - V(tt);
    alpha(tt) = min(alpha0*G(tt), 1);

    V(tt+1) = V(tt) + alpha(tt)*rpe(tt);
end

G = [nan, G];

itiOpt = D./V;
itiMdl = generate_wts(itiOpt, 1, 'logn', 1.5);
end