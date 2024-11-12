function [Xsm] = causal_smooth(X,win)
 %simple code to do smoothing using data from future only
 Xsm = zeros(size(X));

 ns = numel(X);

 for j = 1:ns
     if j < win
         Xsm(j) = 1/j*sum(X(1:j), 'omitnan');    
     else
        Xsm(j) = 1/win*sum(X(j-(win-1):j), 'omitnan');
     end
 end
 