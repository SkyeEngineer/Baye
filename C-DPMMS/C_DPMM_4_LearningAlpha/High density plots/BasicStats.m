function statsmat = BasicStats(data,datasubset,tagnames)
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
% BASICSTATS: Does basic statistical evaluation of a data set. 
% Call it like this:
%
% BasicStats(filename,datasubset,tagnames)
% For example:
% statsmat = BasicStats('SEAsiaPVmat',[1:450]);
% For each tag, the routinr prints the mean value, standard deviation and 
% standard deviation as a percentage of mean. 
%
% Input arguments:
% data: data to be visualized with the time trends arranged as columns. 
% Datasubset: The samples to be examined. Default is all the samples.
% Tagnames: A optional cell array containing tagnames. 
%
% Output argument:
% statsmat: Matrix with m rows where m is the number of tags in the data set. 
% Its columns are: Co1 1 - mean value, Col 2 - standard deviation, 
% Col 3 - standard deviation as percent of mean value. 

%========================================================================

%--------------------------------------------------------------------
% THIS SECTION ORGANISES DATA

% eval(['global ' filename])  % make it available for subsequent calls
%       % NB To clear a global variable type clear global at MATLAB propmt, not just clear.  
% eval(['if isempty(' filename  '), load ' filename '.txt, end'])
% 
% eval(['data = ' filename ';'])

[NR,NC] = size(data);

% Default values for arguments
if (nargin < 2),
   range = [2 inf];  datasubset = 1:NR;
end
if (nargin < 3),
   tagnames = cellstr(num2str([1:NC]'));  % number the tags if no names available
end

data = data(datasubset,:);
[NR,NC] = size(data);

%--------------------------------------------------------------------
% STATS CALCULATIONS
meanvals = mean(data)'; % set point in engineering units
stdvals = std(data)';  % standard deviation of pv in engineering units
stdperc = 100*stdvals./meanvals; % pv standard deviation as percent of mean

%-----------------------
% These lines added 4th Sept 2003. They work out the appropriate number of significant figures for the report
% The std results are rounded to two sig figs and the mean values are rounded to match
% for instance, 14195.53 +/- 1487.42  becomes 14200 +/- 1500. 
stdvalslog = ceil(log10(stdvals));  %find out how many numerals are present in front of the decimal point
stdvals = round(100*stdvals./(10.^stdvalslog)); % numbers in stdvals2sf are now all like **.00
stdvals =  stdvals.*(10.^stdvalslog)/100;

% Round the meanvals to the same precision
meanvals = round(100*meanvals./(10.^stdvalslog));  %find out how many numerals are present in front of the decimal point
meanvals = meanvals.*(10.^stdvalslog)/100;

% Round the percentages to 2 s.f.
stdperclog = ceil(log10(stdperc));  %find out how many numerals are present in front of the decimal point
stdperc = round(100*stdperc./(10.^stdperclog)); % numbers in stdvals2sf are now all like **.00
stdperc =  stdperc.*(10.^stdperclog)/100;
%-----------------------

%--------------------------------------------------------------------
% SHOW RESULTS
statsmat = [meanvals stdvals stdperc];
format short g  % choose a flexible format

disp([' ']) % blank line
disp(['Tag       mean           std           std%'])  
for i = 1:NC,
    disp(tagnames{i})
    disp([statsmat(i,:)])
end
format short  % restore normal format
%==============================================================================

