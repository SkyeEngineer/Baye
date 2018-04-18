function [eta,alpha]=SamplingAlpha(OldAlpha,k,n)
a=0.001;
b=0.001;
eta=betarnd(OldAlpha+1,n);
numerator=a+k-1;
denominator=n*(b-log(eta))-numerator;
pi_eta=numerator/denominator;
gamma1=gamrnd(a+k,b-log(eta));
gamma2=gamrnd(a+k-1,b-log(eta));
alpha=pi_eta*gamma1+(1-pi_eta)*gamma2;


end
