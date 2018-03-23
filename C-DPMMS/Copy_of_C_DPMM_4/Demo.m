clc;clear;close all



addpath('/Users/congtian/Downloads/Project Code/Data');

% % addpath('/Users/congtian/Downloads/Project Code/CVACaseStudy/CVACaseStudy/CVA function')

% % load FaultyCase1
% % load FaultyCase2
% % load FaultyCase3
% % load FaultyCase4
% % load FaultyCase5
% % load Training
load labelleddata
%
vIndex = 1:17; % measurement index
% % X1 = T1(:,vIndex);
% % X2 = T2(:,vIndex);
% % X3 = T3(:,vIndex);
% % XT = Set1_2(:,vIndex);
Data=[Normal_120air_01water' Normal_150air_05water' ];
len=1:1000;%1:length(Data);
da=1;dw=6;% dim index of water and air
ndim=vIndex;

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
Data=Data(ndim,len);
Max=max(Data'); Min=min(Data');
numerator=bsxfun(@minus, Data, Min');
denominator=Max-Min;
y=bsxfun(@rdivide,numerator,denominator');
% % y=Data;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
[row,col]=size(y);


% % % %%  DPMM CLustering
% % % % Parameters of the base distribution G0
% % % % The parameter of Normal Inverse Wishart
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
%% Clustering normal data
% % % [c_st, record, similarity,p_mu,p_sig,p_isig,ind1,cls_size] = gibbsDPM(Data,da,dw,y, hyperG0, alpha, niter, type_algo, doPlot);
% % % % % % Plot air and water flow
% % % c_est0=record(:,1);
% % % c_est=record(:,2);
% % % figure('Name','Air and Water flow plots','NumberTitle','off');
% % % yyaxis left
% % % plot(1:length(Data),Data(da,:));
% % % xlabel('Time(s)')
% % % ylabel('Air Flow')
% % % title('Air and Water flow data')
% % % yyaxis right
% % % plot(1:length(Data),Data(dw,:));
% % % ylabel('Water Flow')
% % % hold on
% % % plot(1:length(c_est),c_est,'k')
% % % ylabel('Water Flow')
% % % legend('Air flow','Water flow','Clusters')
% % % hold off
% % % print('cls','-depsc')
% % % 
% % % save('ind1.mat','ind1');
% % % save('p_mu.mat','p_mu');
% % % save('p_sig.mat','p_sig');
% % % save('p_isig.mat','p_isig');
% % % save('cls_size.mat','cls_size');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% three y axises plot
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








load cls_size.mat
load ind1.mat
load p_isig.mat
load p_mu.mat
load p_sig.mat

FaultName={'Blockage_120air_01water','Blockage_150air_05water',...
   'Diverted-120air_01water','Diverted-150air_05water',...
   'Leakage-120air_01water','Leakage-150air_05water',...
   'NormalSlugging'};
TeX=NormalSlugging(len,ndim)';
XT = NormalSlugging(len,vIndex);

%%%%%%%%%%%%%%Detection%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numerator1=bsxfun(@minus, TeX(ndim,:), Min');
Te=bsxfun(@rdivide,numerator1,denominator');
[dim, ns]=size(Te);
% % Te=TeX(ndim,:);
iter=20;
cl=zeros(ns,iter);
rs=zeros(1,ns);
% %     figure('Name','Fault Detection Results')
    disp('Detecting..........')
for i=1:ns
    for j=1:iter
         cl(i,j) = detection(ind1,alpha, Te(:,i), hyperG0, p_mu, p_sig,p_isig,cls_size);
    end 

    rs(i)=~ismember(mode(cl(i,10:end),2),ind1);
% %     hold on;
% %     plot(i,rs(i),'or','MarkerSize',5,'MarkerFaceColor','r')
% %     xlabel('Time(s)');ylabel('Results')
% %     names = {'Normal'; 'Abnormal';};
% %     set(gca,'ytick',[0,1],'yticklabel',names)
% %     ylim([0-0.5 1+0.5]);
% %     pause(.001) 
end


%% CVA-case study fault detection

%% Section 1 Parameters

cva_alpha = 0.99;  % confidence level
n=2; %n = 25;        % retained state dimension
p = 15;        % length of past observation
f = 15;        % length of future observation
X1 = Normal_120air_01water(:,vIndex);
X2 = Normal_150air_05water(:,vIndex);

%% Section 4 Fault detection result Figure 7(a) and (b) in the paper
[T2mon,Qmon,Ta,Qa]=cvatutor(cva_alpha, n,p,f,X1,X2,XT);

%% CVA PLOT RESULTS %%%%%%%%%%
N=size(T2mon,2);
figure;
subplot(2,1,1),semilogy(1:N,T2mon,'b',[1 N],[Ta Ta],'r-.','linewidth',2); 
ylabel('T^2'); 
xlabel('Sample Number'); 

subplot(2,1,2),semilogy(1:N,Qmon,'b',[1 N],[Qa Qa],'r-.','linewidth',2); 
ylabel('SPE'); 
xlabel('Sample Number'); 
saveas(gcf,'CVA.png')

T1=find(T2mon>=Ta);T2=find(T2mon<Ta);
Q1=find(Qmon>=Qa);Q2=find(Qmon<Qa);
T2mon(T1)=1;T2mon(T2)=0;
Qmon(Q1)=1;Qmon(Q2)=0;
%% C-DPMMs PLOT RESULTS
figure('Name','Fault Detection Results');
subplot(3,1,1)
plot(1:length(rs),rs,'-r','Linewidth',3);
xlabel('Time(s)');ylabel('Results')
names = {'Normal'; 'Abnormal';};
set(gca,'ytick',[0,1],'yticklabel',names)
ylim([0-0.5 1+0.5]);
title('C-DPMMs');

subplot(3,1,2)
plot(1:length(T2mon),T2mon,'-r','Linewidth',3);
xlabel('Time(s)');ylabel('T2')
names = {'Normal'; 'Abnormal';};
set(gca,'ytick',[0,1],'yticklabel',names)
ylim([0-0.5 1+0.5]);
title('CVA-T2')

subplot(3,1,3)
plot(1:length(Qmon),Qmon,'-r','Linewidth',3);
xlabel('Time(s)');ylabel('SEP')
names = {'Normal'; 'Abnormal';};
set(gca,'ytick',[0,1],'yticklabel',names)
ylim([0-0.5 1+0.5]);
title('CVA-Q')

% % print('Fault Detection Results','-depsc')
saveas(gcf,'C-DPMMs.png')

