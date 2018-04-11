clc;clear;close all
% % % clustering FlowConditon From Ruomu
cd ..
addpath(genpath(pwd));
addpath(genpath(pwd),'helpfunctions');
addpath(genpath(pwd),'data');
addpath(genpath(pwd),'High density plots')
load labelleddata.mat
load FlowCondition.mat


load Day2.mat
D2_1=Day2(5161:5761,:); %13:16:00-13:26:00 Normal_120air_01water
D2_2=Day2(10261:10801,:); %14:41:01-14:50:00 Normal_150air_05water


load Day3.mat

D3_1=Day3(1501:2401,:); %10:25:00-10:40:00 Normal_120air_01water

D3_2=Day3(5761:6901,:); %11:36:00-11:55:00 Normal_150air_05water

D3_3=Day3(18301:18601,:);%15:05:00-15:10:00 Normal_150air_05water_1






%% Section 1 Parameters settings
vIndex = [2,3,5:17]; % measurement index
da=1;dw=6;% dim index of water and air



%% Data
%%% Training and Test Data

Training=[FlowCondition(1:end,:)'];

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


%% three y axises plot
figure(3)
ylabels{1}='Air flow';
ylabels{2}='Water flow';
ylabels{3}='Cluster lables';
[ax,hlines] = plotyyy(1:length(Data),Data(8,:),1:length(Data),Data(9,:),1:length(Data),c_est,ylabels)
save('DPMMC_results.mat','c_est');

ncls=length(unique(c_est));
cs=cell(2,ncls);
for i=1:ncls
    cs{1,i}=find(c_est==i);
    cs{2,i}=length(cs{1,i});
end
% % sort(cs{2,:},'descend')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% saving training results

save('ind1.mat','ind1');
save('p_mu.mat','p_mu');
save('p_sig.mat','p_sig');
save('p_isig.mat','p_isig');
save('cls_size.mat','cls_size');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = 'FlowCondition.mat'; % name of the mat file
% % data = importdata(filename); 
tagnames = textread('tagnames.txt','%q'); % Tagnames is a cell array

%--------------------------------------------------------------------
% Specify the range of filter used and samples to be visualized 
[NR,NC] = size(Training');
range = [2 inf];  % no filtering (asks for 2 samples per cycle to d.c.
datasubset = 1:NR; % Plot all data samples

%--------------------------------------------------------------------
% DO THE PLOTS
HighDensityPlot(Training',filename,range,datasubset,tagnames,c_est);

%--------------------------------------------------------------------
% DO THE STATS CALCULATIONS
BasicStats(Training',datasubset,tagnames); 

