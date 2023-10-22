function [tau, newy] = get_tau(xvals, yvals)

%myx = floor(length(xvals)/3);
iters = 25;

if sum(~isnan(yvals))>0
    b = robustfit(xvals(1:3), yvals(1:3));
    %need a tight window for the opt-outs
    if b(2)<0 %if negative slope
        m = min(yvals);
        y = yvals-m;
        [~, bestfit, ~] = fit_exp_decay(xvals, y, iters);
        newy = bestfit(1).*exp(-xvals./bestfit(2)) + m;
    else
        y = -1.*yvals;
        m = min(y);
        y = y - m;
        [~, bestfit, ~] = fit_exp_decay(xvals, y, iters);
        newy = bestfit(1).*exp(-xvals./bestfit(2)) + m;
        newy = newy.*-1;
    end
    tau = bestfit(2);
else
    tau = nan;
    newy = nan(1, length(yvals));
end

end