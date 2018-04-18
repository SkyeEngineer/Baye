function [c_st, R_mu, R_sigma, Ri_sigma] = gibbsDPM_algo2(Data,da,dw, y, hyperG0, alpha, niter, doPlot,learning)

% Gibbs sampler for Dirichlet Process Mixtures of Gaussians
% The base distribution G0 is Normal inverse Wishart of parameters given by
% hyperG0
% Reference: Algorithm 2 of Neal

if doPlot
    figure('name','Gibbs sampling for DPM');
    colormap('default')
    cmap = colormap;
end


[p, n] = size(y);
c_st = zeros(n, niter/2);
U_mu = zeros(p, n);
U_Sigma = zeros(p, p, n);
iU_Sigma = zeros(p, p, n);
ALPHA=[];
ETA=[];
KK=[];
% U_SS is a structure array where U_SS(k) contains the sufficient
% statistics associated to cluster k
U_SS = struct('mu', cell(n, 1), 'kappa', cell(n, 1), ...
    'nu', cell(n, 1), 'lambda', cell(n, 1));

m=zeros(1,n);
c = zeros(n, 1);
% Initialisation
for k=1:n
    c(k) = randi(n); ceil(30*rand); % Sample new allocation uniform
    m(c(k)) = m(c(k)) + 1;
    if m(c(k))>1
        U_SS(c(k)) = update_SS(y(:,k), U_SS(c(k)));
    else
        U_SS(c(k)) = update_SS(y(:,k), hyperG0);
    end
end
ind = unique(c);
for j=1:length(ind)
    [U_mu(:, ind(j)), U_Sigma(:, :, ind(j)), iU_Sigma(:, :, ind(j)) ] = normalinvwishrnd(U_SS(ind(j)));
end

R_mu=[];
R_sigma=[];
Ri_sigma=[];
% Iterations
for i=2:niter
    % Update cluster assignments c
    for k=1:n        
        m(c(k)) = m(c(k)) - 1;
        U_SS(c(k)) = downdate_SS(y(:,k),U_SS(c(k)));
        c(k) = sample_c(c,k,m, alpha, y(:,k), hyperG0, U_mu, U_Sigma,iU_Sigma);
        m(c(k)) = m(c(k)) + 1;
        if m(c(k))>1
            U_SS(c(k)) = update_SS(y(:,k), U_SS(c(k)));
            [U_mu(:, c(k)), U_Sigma(:, :, c(k)), iU_Sigma(:, :, c(k))] = normalinvwishrnd(U_SS(c(k)));
        else
            U_SS(c(k)) = update_SS(y(:,k), hyperG0);
           [U_mu(:, c(k)), U_Sigma(:, :, c(k)), iU_Sigma(:, :, c(k))] = normalinvwishrnd(U_SS(c(k)));
        end       
      
        if doPlot==1
            some_plot(Data,y, U_mu, m, c, k, i, cmap)
        end
    end
    % Update cluster locations U
% %     ind = find(m);
% %     for j=1:length(ind)
% %         [U_mu(:, ind(j)), U_Sigma(:, :, ind(j)), iU_Sigma(:, :, ind(j)) ] = normalinvwishrnd(U_SS(ind(j)));
% %     end
    % Update alpha
% %     nofcluster=length(find(m~=0));
switch(learning)
    case '1'
        [eta, newalpha]=SamplingAlpha(alpha,length(unique(c)),n);
        alpha=newalpha;
        ALPHA=[ALPHA newalpha];
        ETA=[ETA eta];
        KK=[KK length(unique(c))];
    case '0'
        alpha=alpha;
end

    
    fprintf('Iteration %d/%d\n', i, niter)
    fprintf('%d clusters\n', length(unique(c)))
    fprintf('value of alpha %f\n',alpha)
 

    
    if doPlot==2
        some_plot(Data,da,dw,y, U_mu, m, c, k, i, cmap)
    end
    
    if i>niter/2
        c_st(:, i-niter/2) = c;
        R_mu=[R_mu;U_mu];
        R_sigma=[R_sigma; U_Sigma];
        Ri_sigma=[Ri_sigma; iU_Sigma];
    end
    if i==niter
        PROB=[];
        for ii=1:length(ALPHA)
            temp=0;
            for j=1:length(ALPHA)
                probability=MonteCarloAverage(ALPHA(ii),0.001,0.001,ETA(j),n,KK(j));
                temp=temp+probability;
            end
            PROB=[PROB temp/length(ALPHA)];
        end
        figure('Name','Posterior of ALpha')
        [alpha_value,alpha_index]=sort(ALPHA(:));
        
        plot(alpha_value,PROB(alpha_index)/sum(PROB(alpha_index)),'Linewidth',3);
    end
    

end
end


% % function K = sample_c(rcd,idx, m, alpha, z, hyperG0, U_mu, U_Sigma,iU_Sigma)
% % 
% %     c = find(m~=0); % gives indices of non-empty clusters
% %     r = sum(m);
% %     n = m(c).*exp(loggausspdf(repmat(z, 1, length(c))', U_mu(:, c)', U_Sigma(:, :, c),iU_Sigma(:,:,c))');
% % 
% %     n0 = pred(z, hyperG0);
% %     const = sum(n) + alpha*n0;
% % 
% %     p0 = alpha*n0/const; % probability of sampling a new item
% % 
% %     u=rand(1);
% %     if u<p0
% %         K = find(m==0, 1 );
% %     else
% %         u1 = (u-p0);
% %         temp=sort(n,'descend');%% test
% %         ind = find(cumsum(n/const)>=u1, 1 );  
% %         if idx==1
% %             K = c(ind);
% %         else
% %             if randi(2)==1
% %                 K = c(ind);
% %             else
% %                 K = rcd(idx-1);
% %             end
% %         end
% %         
% % 
% %     end
% % end

function some_plot(rawdata,da,dw, z, U_mu, m, c, k, i, cmap)
        ind=find(m);
%         m(ind)
        hold off
%         length(ind)
% % subplot(2,1,1)
        for j=1:length(ind)
            plot(rawdata(da,c==ind(j)),rawdata(dw,c==ind(j)),'.','color',cmap(mod(5*ind(j),63)+1,:), 'markersize', 15);        
            hold on
% %             plot(U_mu(1,ind(j)),U_mu(2,ind(j)),'.','color',cmap(mod(5*ind(j),63)+1,:), 'markersize', 30); 
% %             plot(U_mu(1,ind(j)),U_mu(2,ind(j)),'ok', 'linewidth', 2, 'markersize', 10); 
        end
% %         plot(z(1,k),z(2,k),'or', 'linewidth', 3)            
        title(['i=' num2str(i) ',  k=' num2str(k) ', Nb of clusters: ' num2str(length(ind))]);
        xlabel('Air flow')
        ylabel('Water flow')
        % %             xlim([-3 3]);
        % %             ylim([-3 3]);
        minair=min(rawdata(da,:));minwater=min(rawdata(dw,:));
        maxwair=max(rawdata(da,:));maxwater=max(rawdata(dw,:));
        xlim([minair maxwair]);
        ylim([minwater maxwater]);
        pause(.01)
        
% % hold off       
% % subplot(2,1,2)
% %       
% % % % figure('Name','Air and Water flow plots','NumberTitle','off');
% % yyaxis left
% % plot(1:length(rawdata),rawdata(dw,:));
% % xlabel('Time(s)')
% % ylabel('Water Flow')
% % title('Air and Water flow data')
% % hold on
% % yyaxis right
% % plot(1:length(rawdata),rawdata(da,:));
% % ylabel('Air Flow')
% % hold on
% % plot(1:length(c),c,'k')
% % ylabel('Water Flow')
% % legend('Air flow','Water flow','Clusters')
% % pause(5)
end
