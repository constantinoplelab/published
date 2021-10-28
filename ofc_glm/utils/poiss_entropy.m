function H = poiss_entropy(lambda,kmax)
%calculates entropy of a poisson distribution with rate lambda and max
%number of events k
%
%inputs: lambda: n x 1 set of rate parameters
%        kmax: max number of events that can occur in single time bin

Hsum = 0;
for j = 1:kmax
    Hsum = Hsum + lambda.^(j)*log(factorial(j))/factorial(j);
end

H = lambda.*(1-log(lambda))+ exp(-lambda).*Hsum;

    
    