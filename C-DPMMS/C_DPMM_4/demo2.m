clc;clear;close all
% % % clustering normal data
addpath('/Users/congtian/Downloads/Project Code/Data');
load labelleddata.mat


load Day2.mat
D2_1=Day2(5161:5761,:); %13:16:00-13:26:00 Normal_120air_01water
D2_2=Day2(10261:10801,:); %14:41:01-14:50:00 Normal_150air_05water


load Day3.mat

D3_1=Day3(1501:2401,:); %10:25:00-10:40:00 Normal_120air_01water

D3_2=Day3(5761:6901,:); %11:36:00-11:55:00 Normal_150air_05water

D3_3=Day3(18301:18601,:);%15:05:00-15:10:00 Normal_150air_05water_1




%%
da=1;dw=6;% dim index of water and air
%% Section 1 Parameters settings
vIndex = 1:17; % measurement index
%%% Section 1.1 CVA parameter
cva_alpha = 0.99;  % confidence level
n=2; %n = 25;        % retained state dimension
p = 15;        % length of past observation
f = 15;        % length of future observation


%% Data
%%% Training and Test Data
X1 = Normal_120air_01water(:,vIndex); % operation 1
X2 = Normal_150air_05water(:,vIndex); % operation 2
Training=[D3_2' D2_1' D3_1' D2_2' D3_3'];

FaultName={'Blockage_120air_01water','Blockage_150air_05water',...
   'Diverted_120air_01water','Diverted_150air_05water',...
   'Leakage_120air_01water','Leakage_150air_05water',...
   'NormalSlugging'};
len=1:length(Diverted_120air_01water);
TeX=Diverted_120air_01water(len,vIndex)';
%%% Max-Min normalization
Data=Training(vIndex,:);
Max=max(Data'); Min=min(Data');
numerator=bsxfun(@minus, Data, Min');
denominator=Max-Min;
y=bsxfun(@rdivide,numerator,denominator');
[row,col]=size(y);

%%% Section 1.2 The parameter of Normal Inverse Wishart
% % % % Parameters of the base distribution G0
hyperG0.mu = zeros(row,1);%[0;0];
hyperG0.kappa = 1;
hyperG0.nu = row+2;% hyperG0.nu should be greater than the number of variables set it as 4
hyperG0.lambda = eye(row);

alpha =10;% Scale parameter of DPM
niter = 20; % Number of iterations
doPlot = 2;% do some plots
type_algo = 'CRP'; % other algorithms: 'CRP', 'collapsesCRP', 'slicesampler'

numerator1=bsxfun(@minus, TeX(vIndex,:), Min');
Te=bsxfun(@rdivide,numerator1,denominator');
[dim, ns]=size(Te);


%% C-DPMMs: Clustering based on Dirichlet Process Mixture Mdoels
% % %%% Clustering normal data
[c_st, record, similarity,p_mu,p_sig,p_isig,ind1,cls_size] = gibbsDPM(Training,da,dw,y,hyperG0, alpha, niter, type_algo, doPlot);
 
% % % % % Plot air and water flow
c_est0=record(:,1);
c_est=record(:,2);
figure('Name','Air and Water flow plots','NumberTitle','off');
yyaxis left
plot(1:length(Training),Training(da,:));
xlabel('Time(s)')
ylabel('Air Flow')
title('Air and Water flow data')
yyaxis right
plot(1:length(Training),Training(dw,:));
ylabel('Water Flow')
hold on
plot(1:length(c_est),c_est,'k')
ylabel('Water Flow')
legend('Air flow','Water flow','Clusters')
hold off
print('cls','-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% saving training results
save('ind1.mat','ind1');
save('p_mu.mat','p_mu');
save('p_sig.mat','p_sig');
save('p_isig.mat','p_isig');
save('cls_size.mat','cls_size');
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




%%%%%%%%%%%%%%Detection%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load cls_size.mat
load ind1.mat
load p_isig.mat
load p_mu.mat
load p_sig.mat

cl=zeros(ns,niter);
rs=zeros(1,ns);
% %     figure('Name','Fault Detection Results')
    disp('Detecting..........')
for i=1:ns
    for j=1:niter
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

[T2mon,Qmon,Ta,Qa]=cvatutor(cva_alpha, n,p,f,X1,X2,TeX');

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

%% Detection Rate
CRate=length(find(rs==1))/length(rs);
T2Rate=length(find(T2mon==1))/length(T2mon);
QRate=length(find(Qmon==1))/length(Qmon);

disp(['Detection rates of C-DPMM, T2 and Q are respectively:              '...
    num2str(CRate),'   ', num2str(T2Rate),'   ', num2str(QRate)]);
