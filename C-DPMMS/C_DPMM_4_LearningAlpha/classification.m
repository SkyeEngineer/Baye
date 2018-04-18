function   [K]= classification(norm_cls, alpha, x, hyperG0, U_mu, U_Sigma,iU_Sigma,cls_size,  F_mu, F_Sigma,iF_Sigma,initial,C)

    n = cls_size'.*exp(loggausspdf(repmat(x, 1, length(norm_cls))', U_mu(:, norm_cls)', U_Sigma(:, :, norm_cls),iU_Sigma(:,:,norm_cls))');

    n0 = pred(x, hyperG0);
    const = sum(n) + alpha*n0;

    p0 = alpha*n0/const; % probability of sampling a new item

    u=rand(1);
    if u<p0
        [K]=FaultType(alpha,x,F_mu, F_Sigma,iF_Sigma,initial,C);
       
    else
        u1 = (u-p0);
        ind = find(cumsum(n/const)>=u1, 1 );
        K=-1*norm_cls(ind);

       
    end
end

function [K]=FaultType(alpha,y,F_mu, F_Sigma,iF_Sigma,initial,C)
if length(C)==0
    K=1;
%     F_SS(c(k)) = update_SS(y, initial);
% %        C=[C;0];
% %     C(K)=C(K)+1;
else
    n = C'.*exp(loggausspdf1(repmat(y, 1, length(C))', F_mu(:, 1:length(C))', F_Sigma(:, :, 1:length(C)),iF_Sigma(:,:,1:length(C)))');
    n0 = pred(y, initial);
    const = sum(n) + alpha*n0;
    p0 = alpha*n0/const; % probability of sampling a new item
    
    u=rand(1);
    if u<p0
        K = length(C)+1;
        
    else
        u1 = (u-p0);
        
        ind = find(cumsum(n/const)>=u1, 1 );
        K=ind;
    end
    

end

end