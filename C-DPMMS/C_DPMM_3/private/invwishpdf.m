function pdf=invwishpdf(S,n,lambda)

% pdf d'une distribution inverse Wishart
% Voir définition dans le rapport (différentes définitions possibles)
[p p2]=size(lambda);

if p~=p2
    fprintf('Erreur : Matrice non carrée\n');
end

logtwo=n*p/2*log(2);
logpi=p*(p-1)/4*log(pi);
logdetlambda=-n/2*log(det(lambda));
logdetS=-(n+p+1)/2*log(det(S));
logexptrace=-.5*trace(inv(lambda)*inv(S));

plst=1:p;
gamln=gammaln((n+1-plst)/2);
sumgamln=sum(gamln);

logpdf= logdetlambda + logdetS +logexptrace - (logtwo + logpi + sumgamln);

pdf=exp(logpdf);


             
             
             
             
             
             
             

         
             


