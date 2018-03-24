function [mu, S, iS] = normalinvwishrnd(hyper)

% Sample from a normal inverse Wishart distribution whose parameter are
% given by the structure hyper

mu0 = hyper.mu;
kappa0 = hyper.kappa;
nu0 = hyper.nu;
lambda0 = hyper.lambda;

% Sample S from an inverse Wishart distribution
[iS S] = invwishrnd(nu0,lambda0); 
% Sample mu from a normal distribution
mu = mu0 + chol(S/kappa0)' * randn(length(mu0),1);  
