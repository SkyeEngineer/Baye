function logpdf1 = loggausspdf(x, mu, Sigma, iSigma)

[n, d] = size(x);

[n2,d2] = size(mu);
% Check dimension and center data
if d2 ~= d 
    error(message('X and mu must have the same nb of columns'));
elseif n2 == n 
    x0 = x - mu;
elseif n2 == 1 % mean is a single row, rep it out to match data
    x0 = bsxfun(@minus,x,mu);
else % sizes don't match
    error(message('mu must has 1 row or the same nb as x'));
end



if ndims(Sigma)==2
    [R,err] = chol(Sigma);
    % Standardized data
    xstand =  x0 / R;
    if err~= 0
        error('Sigma must be semi-definite positive');
    end
    % log(sqrt(Sigma))
    logSqrtDetSigma = sum(log(diag(R)));
    
    
    
elseif ndims(Sigma)==3
    xstand = zeros(n,d);
    xstand1 = zeros(n,1);%test
    logSqrtDetSigma = zeros(n,1);
    logSqrtDetSigma1 = zeros(n,1);%test
    for i = 1:n
% %         [R,err] = chol(Sigma(:,:,i));
% %         if err ~= 0
% %             error(message('stats:mvnpdf:BadMatrixSigmaMultiple'));
% % % %         end
% %         xstand(i,:) = x0(i,:) / R;
        xstand1(i,:) = x0(i,:)*iSigma(:,:,i)*x0(i,:)';%test
% %         logSqrtDetSigma(i) = sum(log(diag(R)));
        logSqrtDetSigma1(i)=-log(det(iSigma(:,:,i)))/2;%%test logSqrtDetSigma=logSqrtDetSigma1
    end
end

% Quadratic form

% % xquad = sum(xstand.^2, 2);%xquad=xstand1
% % logpdf = -0.5*xquad - logSqrtDetSigma - d*log(2*pi)/2;
logpdf1= -0.5*xstand1 - logSqrtDetSigma1 - d*log(2*pi)/2;