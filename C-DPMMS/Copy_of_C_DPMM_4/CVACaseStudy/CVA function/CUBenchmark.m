%% A Benchmark Case for Statistic Process Monitoring - Cranfield Multiphase Flow Facility
% This package contains data sets collected in the Cranfield Multiphase
% Flow Facility aiming to serve as a benchmark case for statistic process
% monitoring. Details of the benchmark case are presented in [1]. Some
% examples of Canonical Variate Analysis (CVA) on these data sets are
% presented as follows. 
%% Reference
% [1] C. Ruiz-Cárcel, Y. Cao, D. Mba, L.Lao and R. T. Samuel, Statistical
% Process Monitoring of a Multiphase Flow Facility, _Control Engineering
% Practice_ , V. 42, PP. 74–88, 2015, <http://www.sciencedirect.com/science/article/pii/S0967066115000866#>  

%% Section 1 Parameters
%
clc;close all;clear
vIndex = 1:17; % measurement index
alpha = 0.99;  % confidence level
n=2; %n = 25;        % retained state dimension
p = 15;        % length of past observation
f = 15;        % length of future observation

%% Section 2 Training with normal data sets 2 and 3
% This block loads the available training data sets and selects
% the two data sets used to train the algorithm. The measurements included
% in both data sets are selected according to vIndex
load labelleddata
X1 = Normal_120air_01water(:,vIndex);
X2 = Normal_150air_05water(:,vIndex);

%% Section 3 Fault detection for Faulty Case 3
% This block loads the data set used for the first monitoring example,
% including the same measurements as the training data sets (vIndex)
% % load FaultyCase2
XT = Diverted_150air_05water(:,vIndex);

%% Section 4 Fault detection result Figure 7(a) and (b) in the paper
[T2mon,Qmon,Ta,Qa]=cvatutor(alpha, n,p,f,X1,X2,XT); 

%% Section 5 Fault detection for Faulty Case 5
% This block loads the data set used for the second monitoring example,
% including the same measurements as the training data sets (vIndex)
% % load FaultyCase5
% % XT2 = Set5_2(:,vIndex);

%% Section 6 Fault detection result Figure 9(a) and (b) in the paper
% % cvatutor(alpha, n,p,f,X1,X2,XT2);