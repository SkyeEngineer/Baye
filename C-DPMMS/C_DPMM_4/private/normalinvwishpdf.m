function p=normalinvwishpdf(mu,S,hyper)


mu0=hyper.mu0;
kappa0=hyper.kappa0;
nu0=hyper.nu0;
lambda0=hyper.lambda0;

% pdf d'une distribution normale inverse Wishart

normalpart=gausspdf(mu-mu0,S/kappa0);
invwishpart=invwishpdf(S,nu0,lambda0);

p=normalpart*invwishpart;