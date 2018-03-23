clc;clear;close all
addpath('/Users/congtian/Downloads/Project Code/Data');

load FaultyCase1
load FaultyCase2
load FaultyCase3
load FaultyCase4
load FaultyCase5
load Training
load labelleddata
vIndex = 1:23; % measurement index
X1 = T1(:,vIndex);
X2 = T2(:,vIndex);
X3 = T3(:,vIndex);
XT = Set1_2(:,vIndex);

%% %%%%%%%%%CVA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Section 1 Parameters
% % clc;close all;clear
% % alpha = 0.99;  % confidence level
% % n=25;          % retained state dimension
% % p = 15;        % length of past observation
% % f = 15;        % length of future observation
% % 
% % % Section 4 Fault detection result Figure 7(a) and (b) in the paper
% % [Tr_2D,Te_2D]=cvatutor(alpha, n,p,f,X1,X2,XT);
% % y=Tr_2D(:,1:10000);
%% PCA keep 0.95
% % xtv=[X1' X2'];
% % [u,v]=pca(xtv);
% % d95=find(cumsum(v)./sum(v)<0.95);
% % x=bsxfun(@minus,xtv,mean(xtv,2)); % remove mean
% % y=u(d95,:)*x;
%% Max-Min normalization
Data=[X1'];
len=1:2000;
da=8;dw=9;% dim index of water and air
ndim=1:23;
Data=Data(ndim,len);
% % Max=max(Data'); Min=min(Data');
% % numerator=bsxfun(@minus, Data, Min');
% % denominator=Max-Min;
% % y=bsxfun(@rdivide,numerator,denominator');
y=Data;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
[row,col]=size(y);


%%  DPMM CLustering
% Parameters of the base distribution G0
% The parameter of Normal Inverse Wishart
hyperG0.mu = zeros(row,1);%[0;0];
hyperG0.kappa = 1;
hyperG0.nu = row+2;% hyperG0.nu should be greater than the number of variables set it as 4
hyperG0.lambda = eye(row);
% Scale parameter of DPM
alpha =10;

% Number of iterations
niter = 20;
% do some plots
doPlot = 2;
type_algo = 'CRP'; % other algorithms: 'CRP', 'collapsesCRP', 'slicesampler'

[c_st, record, similarity,p_mu,p_sig,p_isig,ind1,cls_size] = gibbsDPM(Data,da,dw,y, hyperG0, alpha, niter, type_algo, doPlot);
% % % Plot air and water flow
c_est0=record(:,1);
c_est=record(:,2);
figure('Name','Air and Water flow plots','NumberTitle','off');
yyaxis left
plot(1:length(Data),Data(da,:));
xlabel('Time(s)')
ylabel('Air Flow')
title('Air and Water flow data')
yyaxis right
plot(1:length(Data),Data(dw,:));
ylabel('Water Flow')
hold on
plot(1:length(c_est),c_est,'k')
ylabel('Water Flow')
legend('Air flow','Water flow','Clusters')
hold off
print('cls','-depsc')

% % ylabels{1}='Air flow';
% % ylabels{2}='Water flow';
% % ylabels{3}='Cluster lables';
% % [ax,hlines] = plotyyy(1:length(Data),Data(8,:),1:length(Data),Data(9,:),1:length(Data),c_est,ylabels)
% % save('DPMMC_results.mat','c_est');

% % ncls=length(unique(c_est));
% % cs=cell(2,ncls);
% % for i=1:ncls
% %     cs{1,i}=find(c_est==i);
% %     cs{2,i}=length(cs{1,i});
% % end
% % sort(cs{2,:},'descend')













% % %%%%%%%%%%%%%%Detection%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TeX=XT';
[dim, ns]=size(TeX);
% % numerator1=bsxfun(@minus, TeX(ndim,:), Min');
% % Te=bsxfun(@rdivide,numerator1,denominator');
Te=TeX(ndim,:);
iter=20;
cl=zeros(ns,iter);
rs=zeros(1,ns);
% %     figure('Name','Fault Detection Results')
    disp('Detecting..........')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%randomly select out # samples from each cluster
samples=cell(1,length(ind1));
for i=1:length(ind1)
    m0=find(c_est0==ind1(i));
    n0=length(m0);
    if n0<niter
        % padarray to length=20
        n2=randi(n0,[1 niter-n0]);
        m2=[m0 m0(n2)];
        samples{i}=y(:,m2);
    else
        n1=randi(n0,[1 niter]);
        m1=m0(n1);
        samples{i}=y(:,m1);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:ns
    cl(i,:) = detection(niter,samples,ind1,alpha, Te(:,i), hyperG0, p_mu, p_sig,p_isig,cls_size);
    rs(i)= sum(cl(i,:));
    
% %     rs(i)=~ismember(mode(cl(i,1:10),2),[1:size(ind1,3)]);
% %     hold on;
% %     plot(i,rs(i),'or','MarkerSize',5,'MarkerFaceColor','r')
% %     xlabel('Time(s)');ylabel('Results')
% %     names = {'Normal'; 'Abnormal';};
% %     set(gca,'ytick',[0,1],'yticklabel',names)
% %     ylim([0-0.5 1+0.5]);
% %     pause(.001) 
end
    figure('Name','Fault Detection Results');
    plot(1:length(rs),rs);
    xlabel('Time(s)');ylabel('Results')
    names = {'Normal'; 'Abnormal';};
    set(gca,'ytick',[0,1],'yticklabel',names)
    ylim([0-0.5 1+0.5]);

print('Fault Detection Results','-depsc')