function    R= detection(niter,ss,norm_cls, alpha, x, hyperG0, U_mu, U_Sigma,iU_Sigma,cls_size);

% %     norm_cls = find(m~=0); %  indices of clusters containing normal operations
% %     r = sum(m);
% %     n = cls_size'.*exp(loggausspdf(repmat(x, 1, length(norm_cls))', U_mu(:, norm_cls)', U_Sigma(:, :, norm_cls),iU_Sigma(:,:,norm_cls))');
[~,q]=size(ss);
n=q+1;
[p,~]=size(ss{1});
U_SS = struct('mu', cell(n, 1), 'kappa', cell(n, 1), ...
    'nu', cell(n, 1), 'lambda', cell(n, 1));
U_mu = zeros(p,n);
U_Sigma = zeros(p, p,n);
iU_Sigma = zeros(p, p,n);
R=[];
for i=1:niter
    for j=1:q
        U_SS(j)= update_SS(ss{j}(:,i),hyperG0);
        [U_mu(:, j), U_Sigma(:,:, j), iU_Sigma(:,:, j)] = normalinvwishrnd(U_SS(j));
    end
    U_SS(n)= update_SS(x,hyperG0);
    [U_mu(:, n), U_Sigma(:,:, n), iU_Sigma(:,:, n)] = normalinvwishrnd(U_SS(n));
    pdf = exp(loggausspdf(repmat(x, 1, q)', U_mu(:, 1:q)', U_Sigma(:, :, 1:q),iU_Sigma(:,:,1:q))');
    n0 = pred(x, hyperG0);
    const = sum(pdf) + alpha*n0;
    
    p0 = alpha*n0/const; % probability of sampling a new item
    
    u=rand(1);
    if u<p0
        K=1;
    else
        u1 = (u-p0);
        ind = find(cumsum(pdf/const)>=u1, 1 );
        K=0;
    end
    R=[R K];
end





end