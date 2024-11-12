% Generate Model WTs by drawing from gamma or logn distribution
function wt_mdl = generate_wts(wt_opt, nsims, noise_arg, params)

if strcmpi(noise_arg, 'gamma')
    Sigma = params(1);
    
    k = wt_opt.^2./Sigma;
    theta = Sigma./wt_opt;

    wt_mdl = mean(cell2mat(arrayfun(@(x) gamrnd(k, theta), 1:nsims,...
        'UniformOutput', false)), 2);
elseif strcmpi(noise_arg, 'logn') || strcmpi(noise_arg, 'logn_secondhalf') 
    Sigma = params(1);
    
    var = log(Sigma./wt_opt.^2 + 1);
    mu = log(wt_opt) - (var/2);
    
    wt_mdl = mean(cell2mat(arrayfun(@(x) lognrnd(mu, sqrt(var)),...
        1:nsims, 'UniformOutput', false)), 2, 'omitnan');
elseif strcmpi(noise_arg, 'exg_wt')
    
    sigma = params(1);
    tau = params(2);
    m = wt_opt - tau;
    
    wt_mdl = arrayfun(@(mu) mean(exgrnd(mu, sigma, tau, [1 nsims])), m);
elseif strcmpi(noise_arg, 'exg_iti')
   
    sigma = params(1);
    tau = params(2);
    m = log(wt_opt) - tau;
    
    wt_mdl = arrayfun(@(mu) mean(exgrnd(mu, sigma, tau, [1 nsims])), m);
    wt_mdl = exp(wt_mdl);
elseif strcmpi(noise_arg, 'normal')
    
    Sigma = params(1);
    
    wt_mdl = arrayfun(@(mu) mean(normrnd(mu, Sigma, [1 nsims])), wt_opt);
end

end
