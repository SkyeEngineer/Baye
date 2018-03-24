function    K= detection(norm_cls, alpha, x, hyperG0, U_mu, U_Sigma,iU_Sigma,cls_size);

% %     norm_cls = find(m~=0); %  indices of clusters containing normal operations
% %     r = sum(m);
    n = cls_size'.*exp(loggausspdf(repmat(x, 1, length(norm_cls))', U_mu(:, norm_cls)', U_Sigma(:, :, norm_cls),iU_Sigma(:,:,norm_cls))');

    n0 = pred(x, hyperG0);
    const = sum(n) + alpha*n0;

    p0 = alpha*n0/const; % probability of sampling a new item

    u=rand(1);
    if u<p0
        ab=setdiff(1:length(U_mu),norm_cls);
        K=ab(1);
% %         K = find(m==0, 1 );
    else
        u1 = (u-p0);
% %         temp=sort(n,'descend');%% test
        ind = find(cumsum(n/const)>=u1, 1 );
        K=norm_cls(ind);
% %         if idx==1
% %             K = c(ind);
% %         else
% %             if randi(2)==1
% %                 K = c(ind);
% %             else
% %                 K = c(ind);%rcd(idx-1);
% %             end
% %         end
        

    end
end