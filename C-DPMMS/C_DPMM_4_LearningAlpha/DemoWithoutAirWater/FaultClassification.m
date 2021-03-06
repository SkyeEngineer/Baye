clc;clear;close all
% % % training data combined from day2 & day3 tesing on all types of faults
% % % Variable 3 and 7 are air temperature and water temperature
cd ..
addpath(genpath(pwd));
addpath(genpath(pwd),'helpfunctions');
addpath(genpath(pwd),'data');
addpath(genpath(pwd),'High density plots');
% cd('C:\Users\Tian Cong\Desktop\Github\Bayes\Baye\C-DPMMS\C_DPMM_4\DemoWithoutAirWater')

load labelleddata.mat
load Leakage_120air_01water.mat



%% Section 1 Parameters settings
% vIndex = 1:17;
vIndex = [2,4,5,8:17];% measurement index
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
Fault=[NormalSlugging];
len=1:length(Fault);
TeX=Fault(len,vIndex)';
% % save('TeX.mat','TeX');
%%% Max-Min normalizatio
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

numerator1=bsxfun(@minus, TeX(:,:), Min');
Te=bsxfun(@rdivide,numerator1,denominator');
[dim, ns]=size(Te);


% %% C-DPMMs: Clustering based on Dirichlet Process Mixture Mdoels
% % % % %%% Clustering normal data
% % [c_st, record, similarity,p_mu,p_sig,p_isig,ind1,cls_size] = gibbsDPM(Training,da,dw,y,hyperG0, alpha, niter, type_algo, doPlot);
% % % %  
% % % % % % % % % Plot air and water flow
% % c_est0=record(:,1);
% % c_est=record(:,2);
% % figure('Name','Air and Water flow plots','NumberTitle','off');
% % yyaxis left
% % plot(1:length(Training),Training(da,:));
% % xlabel('Time(s)')
% % ylabel('Air Flow')
% % title('Air and Water flow data')
% % yyaxis right
% % plot(1:length(Training),Training(dw,:));
% % ylabel('Water Flow')
% % hold on
% % plot(1:length(c_est),c_est,'k')
% % legend('Air flow','Water flow','Clusters')
% % hold off
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% saving training results
% % % 
% % save('ind1.mat','ind1');
% % save('p_mu.mat','p_mu');
% % save('p_sig.mat','p_sig');
% % save('p_isig.mat','p_isig');
% % save('cls_size.mat','cls_size');








%%%%%%%%%%%%%%Detection%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load cls_size.mat
load ind1.mat
load p_isig.mat
load p_mu.mat
load p_sig.mat

cl=zeros(ns,niter);
rs=zeros(1,ns);
operation=zeros(1,ns);
F_mu = [];
F_Sigma = [];
iF_Sigma = [];
C=[];
initial.mu = zeros(row,1);%[0;0];
initial.kappa = 1;
initial.nu = row+2;% hyperG0.nu should be greater than the number of variables set it as 4
initial.lambda = eye(row);
F_SS = struct('mu', cell(1, 1), 'kappa', cell(1, 1), ...
    'nu', cell(1, 1), 'lambda', cell(1, 1));
% %     figure('Name','Fault Detection Results')
disp('Detecting..........')
for i=1:ns
    for j=1:niter
         [cl(i,j)] = classification(ind1,alpha, Te(:,i), hyperG0, p_mu, p_sig,p_isig,cls_size,F_mu, F_Sigma,iF_Sigma,initial,C);
    end 
    operation(i)=mode(cl(i,10:end),2);
    rs(i)=~ismember(operation(i),-1*(ind1));% multiplying -1 aims to distuiguish the same index between fault and normal
    if rs(i)==1
        K=operation(i);
        if length(C)<K
            C=[C;0];
            C(K)=C(K)+1;
        else
            C(K)=C(K)+1;
        end
        if C(K)>1
            F_SS(K) = update_SS(Te(:,i), F_SS(K));
        else
            F_SS(K) = update_SS(Te(:,i), initial);
        end
         [F_mu(:, K), F_Sigma(:, :, K), iF_Sigma(:, :, K) ] = normalinvwishrnd(F_SS(K));
    end
    
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
subplot(4,1,1)
plot(1:length(rs),rs,'-r','Linewidth',3);
xlabel('Time(s)');ylabel('Results')
names = {'Normal'; 'Abnormal';};
set(gca,'ytick',[0,1],'yticklabel',names)
ylim([0-0.5 1+0.5]);
title('C-DPMMs');


subplot(4,1,2)
temp=find(operation<0);
operation(temp)=0;
plot(1:length(operation),operation,'-r','Linewidth',3);
xlabel('Time(s)');ylabel('Results')
% names = {'Normal'; 'Abnormal';};
% set(gca,'ytick',[0,1],'yticklabel',names)
ylim([0-0.5 1+0.5]);
title('C-DPMMs');

subplot(4,1,3)
plot(1:length(T2mon),T2mon,'-r','Linewidth',3);
xlabel('Time(s)');ylabel('T2')
names = {'Normal'; 'Abnormal';};
set(gca,'ytick',[0,1],'yticklabel',names)
ylim([0-0.5 1+0.5]);
title('CVA-T2')

subplot(4,1,4)
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




% % %% High Density Plots
% % filename = 'TeX.mat'; % name of the mat file
% % % % data = importdata(filename); 
% % tagnames = textread('tagnote.txt','%q'); % Tagnames is a cell array
% % 
% % %--------------------------------------------------------------------
% % % Specify the range of filter used and samples to be visualized 
% % [NR,NC] = size(TeX');
% % range = [2 inf];  % no filtering (asks for 2 samples per cycle to d.c.
% % datasubset = 1:NR; % Plot all data samples
% % 
% % %--------------------------------------------------------------------
% % % DO THE PLOTS
% % HighDensityPlot(TeX',filename,range,datasubset,tagnames,rs);
% % 
% % %--------------------------------------------------------------------
% % % DO THE STATS CALCULATIONS
% % BasicStats(TeX',datasubset,tagnote); 