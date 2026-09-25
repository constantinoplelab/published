function [Xsm] = causal_smooth_SS(X,win)
 %simple code to do smoothing using data from future only
 Xsm = zeros(size(X));

 ns = numel(X);

 for j = 1:ns
     if isnan(X(j)) 
         Xsm(j) = nan;
     elseif j < win
         Xsm(j) = 1/j*sum(X(1:j));    
     else
        Xsm(j) = 1/win*sum(X(j-(win-1):j));
     end
 end
 