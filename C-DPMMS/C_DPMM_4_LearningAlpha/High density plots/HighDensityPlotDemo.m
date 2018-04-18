function [ dummy] = HighDensityPlotDemo
%========================================================================
% Author:  Nina Thornhill
%          Department of Electronic and Electrical Engineering
%          University College London
%          Torrington Place
%          London WC1E 7JE
%
% revised by Ruomu Tan 20180329
%
%========================================================================
% Date:    September 2001
%========================================================================
% HIGHDENSITYPLOTDEMO loads and plots the time series data obtained in the 
% PRONTO Cranfield data set (FlowCondition.mat)
% The plots produced are time trends, autocovariances and power spectra.
% It also performs basic statistical calculations on the data. 
%========================================================================

%--------------------------------------------------------------------
% Import the data and tag names
filename = 'FlowCondition.mat'; % name of the mat file
data = importdata(filename); 
tagnames = textread('tagnames.txt','%q'); % Tagnames is a cell array

%--------------------------------------------------------------------
% Specify the range of filter used and samples to be visualized 
[NR,NC] = size(data);
range = [2 inf];  % no filtering (asks for 2 samples per cycle to d.c.
datasubset = 1:NR; % Plot all data samples

%--------------------------------------------------------------------
% DO THE PLOTS
HighDensityPlot(data,filename,range,datasubset,tagnames);

%--------------------------------------------------------------------
% DO THE STATS CALCULATIONS
BasicStats(data,datasubset,tagnames); 

end