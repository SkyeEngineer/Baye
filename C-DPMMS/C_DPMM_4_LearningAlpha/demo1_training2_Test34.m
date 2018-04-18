clc;clear;close all
% % % training data from day2 tesing on day3 and day4
% % addpath('/Users/congtian/Downloads/Project Code/Data');
% % load labelleddata.mat
% % load T1.mat
% % load T2.mat
% % load X.mat
addpath(genpath(pwd));
addpath(genpath(pwd),'helpfunctions');
addpath(genpath(pwd),'data');

load Day2.mat
Normal_120air_01water=Day2(5161:5761,:); %13:16:00-13:26:00
Normal_150air_05water=Day2(10261:10801,:); %14:41:01-14:50:00

D2_Blockage_120air_01water=Day2(5762:9841,:); %13:26:01-14:34:00

D2_Blockage_150air_05water=Day2(10802:14581,:); %14:50:01-15:53:00
D2_Leakage_150air_05water=Day2(15001:17641,:); %16:00:00-16:44:00
load Day3.mat

D3_Normal_120air_01water=Day3(1501:2401,:); %10:25:00-10:40:00
D3_Leakage_120air_01water=Day3(2402:5161,:); %10:40:01-11:26:00
D3_Normal_150air_05water=Day3(5761:6901,:); %11:36:00-11:55:00
D3_Leakage_150air_05water=Day3(6902:9421,:);%11:55:01-12:37:00
D3_Diverted_120air_01water=[Day3(10204:13425,:);Day3(17101:18061,:)];%12:50-13:43:00, 14:45:00-15:01:00
D3_Normal_150air_05water_1=Day3(18301:18601,:);%15:05:00-15:10:00
D3_Diverted_150air_05water=Day3(18602: 21901,:); % 15:10:01-16:05:00

Te_day3=[D3_Normal_120air_01water;D3_Leakage_120air_01water;...
         D3_Normal_150air_05water;D3_Leakage_150air_05water; ...
         D3_Diverted_120air_01water; D3_Normal_150air_05water_1 ;...
         D3_Diverted_150air_05water];
load Day4.mat

NormalSlugging=Day4;
  


%% Section 1 Parameters settings
vIndex = 1:17; % measurement index
da=1;dw=6;% dim index of water and air

%%% Section 1.1 CVA parameter
cva_alpha = 0.99;  % confidence level
n=25;        % retained state dimension
p = 15;        % length of past observation
f = 15;        % length of future observation


%% Data
%%% Training and Test Data
X1 = Normal_120air_01water(:,vIndex); % operation 1
X2 = Normal_150air_05water(:,vIndex); % operation 2
Training=[Normal_120air_01water' Normal_150air_05water' ];

FaultName={'Blockage_120air_01water','Blockage_150air_05water',...
   'Diverted_120air_01water','Diverted_150air_05water',...
   'Leakage_120air_01water','Leakage_150air_05water',...
   'NormalSlugging'};
len=1:length(Te_day3);
TeX=Te_day3(len,vIndex)';
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
T2Rate=length(find(T2mon==1))/length(rs);
QRate=length(find(Qmon==1))/length(rs);

disp(['Detection rates of C-DPMM, T2 and Q are respectively:              '...
    num2str(CRate),'   ', num2str(T2Rate),'   ', num2str(QRate)]);
