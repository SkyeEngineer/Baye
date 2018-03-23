function c_st = gibbsDPM_algo4(y, hyperG0, alpha, niter, doPlot)

% Gibbs sampler for Dirichlet Process Mixtures of Gaussians
% The base distribution G0 is Normal inverse Wishart of parameters given by
% hyperG0
% Slice sampler
% Reference: Kalli, Griffin, Walker (2011) and Fall and Barat (2012)

if doPlot
    figure('name','Gibbs sampling for DPM');
    colormap('default')
    cmap = colormap;
end


[p, n] = size(y);
c_st = zeros(n, niter/2);
U_mu = zeros(p, n);
U_Sigma = zeros(p, p, n);

% U_SS is a structure array where U_SS(k) contains the sufficient
% statistics associated to cluster k
U_SS = struct('mu', cell(n, 1), 'kappa', cell(n, 1), ...
    'nu', cell(n, 1), 'lambda', cell(n, 1));
for k=1:n
    U_SS(k) = hyperG0;
end

m=zeros(1,n);
c = zeros(n, 1);
% Initialisation: everybody assigned to a different cluster
for k=1:n
    c(k) = ceil(50*rand); % Sample new allocation uniform in 1:10
    U_SS(c(k)) = update_SS(y(:,k), U_SS(c(k)));
    [U_mu(:, c(k)), U_Sigma(:, :, c(k))] = normalinvwishrnd(U_SS(c(k)));
    m(c(k)) = m(c(k))+1;
end
c_st(:, 1) = c;

% Iterations
for i=2:niter    
    % Update cluster assignments c with slice sampling
    [c, m, U_mu, U_sigma, U_SS] = slicesample_c(y, c, U_mu, U_Sigma, U_SS, hyperG0, alpha);

    % Update cluster locations U
    ind = find(m);
    for j=1:length(ind)
        [U_mu(:, ind(j)), U_Sigma(:, :, ind(j))] = normalinvwishrnd(U_SS(ind(j)));
    end
    
    fprintf('Iteration %d/%d\n', i, niter)
    fprintf('%d clusters\n\n', length(ind))
    
    if i>niter/2
        c_st(:, i-niter/2) = c;
    end

    if doPlot
        some_plot(y, U_mu, m, c, 1, i, cmap)
    end
end

end


function some_plot(z, U_mu, m, c, k, i, cmap)
        ind=find(m);
        m(ind)
        hold off
        length(ind)
        for j=1:length(ind)
            plot(z(1,c==ind(j)),z(2,c==ind(j)),'.','color',cmap(mod(5*ind(j),63)+1,:), 'markersize', 15);        
            hold on
            plot(U_mu(1,ind(j)),U_mu(2,ind(j)),'.','color',cmap(mod(5*ind(j),63)+1,:), 'markersize', 30); 
            plot(U_mu(1,ind(j)),U_mu(2,ind(j)),'ok', 'linewidth', 2, 'markersize', 10); 
        end
        plot(z(1,k),z(2,k),'or', 'linewidth', 3)            
        title(['i=' num2str(i) ',  k=' num2str(k) ', Nb of clusters: ' num2str(length(ind))]);
        xlabel('X')
        ylabel('Y') 
            xlim([-3 3]);
            ylim([-3 3]);
        pause(.01)
end
