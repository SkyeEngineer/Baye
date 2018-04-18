function probability=MonteCarloAverage(alpha,a,b,eta,n,k)
power1=-alpha*(b-log(eta));
power2=a+k-1;
probability=(alpha^(power2))*(exp(power1))+n*(alpha^(power2-1))*(exp(power1));

end