function c_st = gibbsDPM_algo1(y, hyperG0, alpha, niter, doPlot)

% Gibbs sampler for Dirichlet Process Mixtures of Gaussians
% The base distribution G0 is Normal inverse Wishart of parameters given by
% hyperG0
% Reference: Algorithm 1 of Neal

if doPlot
    figure('name','Gibbs sampling for DPM');
    colormap('default')
    cmap = colormap;
end

[p, n] = size(y);
c_st = zeros(n, niter/2);
theta_mu = zeros(p, n);
theta_Sigma = zeros(p, p, n);

% Initialisation
hyper = update_SS(y(:, 1), hyperG0);
[theta_mu(:, 1), theta_Sigma(:,:, 1)] = normalinvwishrnd(hyper); 
for k=2:n
    [theta_mu(:, k), theta_Sigma(:, :, k)] = sample_theta(alpha,...
            y(:,k), hyperG0, theta_mu(:, 1:k-1), theta_Sigma(:, :, 1:k-1));
end

% Iterations
for i=2:niter
    for k=1:n
        % Sample theta_k | theta_{-k}, y_k
        ind_notk = (1:n)~=k;
        [theta_mu(:, k), theta_Sigma(:, :, k)] = sample_theta(alpha,...
            y(:,k), hyperG0, theta_mu(:, ind_notk), theta_Sigma(:, :, ind_notk));
        if doPlot==1
            some_plot(y, theta_mu, k, i, n, cmap);
        end        
    end
    % get the partition
    [~, ~, c] = unique(theta_mu(1, :));

    fprintf('Iteration %d/%d\n', i, niter);
    fprintf('%d clusters\n\n', length(unique(c)));   
    
    if doPlot==2
            some_plot(y, theta_mu, k, i, n, cmap);
    end   

    if i>niter/2
        c_st(:, i-niter/2) = c;
    end
end

end

%%%% Subfunctions

function [theta_mu_k, theta_Sigma_k] = sample_theta(alpha,...
            z, hyperG0, theta_mu_notk, theta_Sigma_notk)
        
    p = size(theta_mu_notk, 2);

    n = exp(loggausspdf(repmat(z, 1, p)', theta_mu_notk', theta_Sigma_notk));
    n0 = pred(z, hyperG0);
    const = alpha*n0 +sum(n);
    p0 = alpha*n0 / const;

    u =rand;
    if u<p0
        % Sample new value
        hyper = update_SS(z, hyperG0);
        [theta_mu_k, theta_Sigma_k] = normalinvwishrnd(hyper);    
    else
        % Sample old value
        u1 = u - p0;
        ind = find(cumsum(n/const)>=u1, 1);
        theta_mu_k = theta_mu_notk(:, ind);
        theta_Sigma_k = theta_Sigma_notk(:, :, ind);   
    end
end



function some_plot(z, theta_mu, k, i, n, cmap)

    hold off
    ind = unique(theta_mu(1, :));
    c = zeros(n, 1);
    U_mu = zeros(2, length(ind));
    for j=1:length(ind)
        ind2 = find(theta_mu(1, :)==ind(j));
        c(ind2) = j;
        U_mu(:, j) = theta_mu(:, ind2(1));
    end
    for j=1:length(ind)
        plot(z(1,c==j),z(2,c==j),'.','color',cmap(mod(5*j,63)+1,:), 'markersize', 15);        
        hold on
        plot(U_mu(1,j),U_mu(2,j),'.','color',cmap(mod(5*j,63)+1,:), 'markersize', 30); 
        plot(U_mu(1,j),U_mu(2,j),'ok', 'linewidth', 2, 'markersize', 10); 
    end
    plot(z(1,k),z(2,k),'or', 'linewidth', 3)            
    title(['i=' num2str(i) ',  k=' num2str(k) ', Nb of clusters: ' num2str(length(ind))]);
    xlabel('X')
    ylabel('Y') 
    title(['i=' num2str(i) ',  k=' num2str(k) ', Nb of clusters: ' num2str(length(ind))]);
    xlabel('X')
    ylabel('Y')  
    xlim([-3 3]);
    ylim([-3 3]);
    pause(.01)
end