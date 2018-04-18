function K = sample_c(rcd,idx, m, alpha, z, hyperG0, U_mu, U_Sigma,iU_Sigma)

    c = find(m~=0); % gives indices of non-empty clusters
    r = sum(m);
    n = m(c).*exp(loggausspdf1(repmat(z, 1, length(c))', U_mu(:, c)', U_Sigma(:, :, c),iU_Sigma(:,:,c))');

    n0 = pred(z, hyperG0);
    const = sum(n) + alpha*n0;

    p0 = alpha*n0/const; % probability of sampling a new item

    u=rand(1);
    if u<p0
        K = find(m==0, 1 );
    else
        u1 = (u-p0);
        temp=sort(n,'descend');%% test
        ind = find(cumsum(n/const)>=u1, 1 );  
        if idx==1
            K = c(ind);
        else
            if randi(2)==1
                K = c(ind);
            else
                K = c(ind);%rcd(idx-1);
            end
        end
        

    end
end