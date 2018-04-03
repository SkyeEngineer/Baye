function [dummy] = HighDensityPlot(data,filename,range,datasubset,tagnames,clusters)

addpath(genpath(pwd));

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
% HIGHDENSITYPLOT: Plots time trends, autocovariances and spectra of data, 
% or a subset of it. The data can be filtered if desired.
% Call it like this:
%
% HighDensityPlot(filename,range,datasubset,tagnames);
% For example:
% HighDensityPlot('SEAsiaPVmat',[2 inf],[100:400]);
%
% Input arguments:
% data: data set to be visualized with time trends arranged as columns.  
% Filename: Filename.txt is the name of the file holding the data with 
% Range: Range of oscillation periods to be detected. For instance [2 1500]
% detectes oscillations having 2 samples per cycle up to 1500 samples per cycle. 
% Default value is [2 inf] if the user does not specify a range.
% Datasubset: The samples to be examined. Default is all the samples.
% Tagnames: A optional cell array containing tagnames. 



%--------------------------------------------------------------------
% THIS SECTION ORGANISES DATA

[NR,NC] = size(data);

% Default values for arguments
if (nargin < 2),
   range = [2 inf];  datasubset = 1:NR;
end
if (nargin < 3),
   datasubset = 1:NR;
end
if (nargin < 4),
   tagnames = cellstr(num2str([1:NC]'));  % number the tags if no names available
end

% datasubset must be even in length
index = [1:2*floor(length(datasubset)/2)]; %guaranteed to have an even number of elements
datasubset = datasubset(index);
data = data(datasubset,:);
[NR,NC] = size(data);

%--------------------------------------------------------------------
% FILTER DESIGN
lenspec = NR;
cutfreq = [1/range(2) min(1/range(1),0.5)];
cutfreq = floor(cutfreq*lenspec)+1; % convert fractions of the
          % sampling frequency to indexes in the spectrum vector. 
          % e.g. if the spectrum has 4096 frequnecy channels zero frequency 
          % maps to 1 and 0.5 maps to the Nyquist frequency at 2049.
cutfreq(1) = cutfreq(1)+not(cutfreq(1)-1); % add 1 if cutfreq(1)=1
% If user specified 0 (d.c.) in cutfreqsUser alter the lowest channel index in
% cutfreq to 2. This change makes it easier to calculate the aliases
% because the d.c. is special, it has no alised counterpart.
% (since d.c. is always removed anyway it makes no difference to the result.)     
% NumRanges = length(cutfreq)-1; % N lampposts means N-1 spaces

% CHECK FILTER SIZES
deltaf = diff(cutfreq);
fo = mean(cutfreq);
% filter widths as percent of the center frequency - should be more than 0.4
filtersize = deltaf./fo;
if filtersize<0.4,
   disp(['Range is too narrow'])
end

%--------------------------------------------------------------------
% CALCULATE FILTERED SPECTRA, FILTERED TRENDS AND FILTERED AUTOCOVARIANCE
maxautocorr = min(floor(max(range)*4),lenspec/2); % Make sure there are at least 5 cycles
                      % and use floor function to ensure maxautocorr is integer.
numfreq = lenspec/2;  % the part of the spectrum for plotting
 
autocorrmat = []; powspecmat = []; trendmat = []; plotspecmat = [];
for i = 1:NC,
   datain = data(:,i);
   % mean centre 
   datain = datain - mean(datain);
   spec = fft(datain);
   spec=spec(:); % force it to be a column vector
   cutfreqlo = cutfreq(1); 
   cutfreqhi = cutfreq(2);
   % Channels to be retained are [cutfreqlo:cutfreqhi-1] and aliased channels 
   % [lenspec+2-(cutfreqhi-1):lenspec+2-cutfreqlo]. 
   % e.g. if cutfreqlo = 2 and cutfreqhi = 100 and the length of the spectrum is
   % 4096 we should retain channels 2:99 and their aliases in 3999:4096.
   % If cutfreqlo = 100 and cutfreqhi = 200 then retain channels 100 to 199
   % and aliased channels 3899:3998.
   specFilt = zeros(1,lenspec); 
   specFilt(cutfreqlo:cutfreqhi-1)=spec(cutfreqlo:cutfreqhi-1);
   specFilt(lenspec+3-cutfreqhi:lenspec+2-cutfreqlo) = ...
                             spec(lenspec+3-cutfreqhi:lenspec+2-cutfreqlo);
%  ----------------------------------------------------------
   % Next two lines added on 6th Sept 2002
   if (range(1) == 2),
       specFilt(cutfreqhi)=spec(cutfreqhi);  % deal with the Nyquist frequency
   end
%  ----------------------------------------------------------
   trend = real(ifft(specFilt)); trend = trend(:);
   trendmat = [trendmat trend];
   powerspec = specFilt.*conj(specFilt);
   % For plotting, scale to the same maximum value
   plotspec = powerspec(1:numfreq)/max(powerspec(1:numfreq)); plotspec = plotspec(:);
   plotspecmat = [plotspecmat plotspec];
   % For analysis, must be normalised
   powerspec = powerspec/sum(powerspec); % normalise
   powerspec = powerspec(:); % force a column vector
   powspecmat = [powspecmat powerspec];
   autocorr = NR*real(ifft(powerspec)); % Have to times by NR to get ACF(0)=1
   autocorr = autocorr(:); % force a column vector
   autocorrmat = [autocorrmat autocorr(1:maxautocorr)];   
end

%--------------------------------------------------------
% PLOT FILTERED TRENDS
% Mean centre and normalise
for i = 1:NC
   warning off  % suppress warnings because some trends have stdev of 0
   trendmat(:,i)=(trendmat(:,i)-mean(trendmat(:,i)))...
      /std(trendmat(:,i));
   warning backtrace  % warnings back on with line number references
end

figure, clf reset, hold on
offset = 5;
axmin = datasubset(1)-1; axmax = datasubset(NR);
% % yyaxis left
for i = 1:NC
   plot(datasubset,trendmat(:,i)-offset*i,'k','linewidth',1)
   text(axmax*1.01,-offset*i,num2str(i),'fontsize',20,'fontname',...
		'helvetica','verticalalignment','middle','horizontalalignment','left')
end
axis([axmin axmax -(NC+1)*offset 0]),grid
% % set(gca,'YTick',[-NC*offset:offset:-offset])
set(gca,'Yticklabel',tagnames(NC:-1:1))
set(gca,'fontsize',14,'fontname','helvetica','linewidth',1.5,'box','on')
xlabel('time/sample interval','fontsize',20,'fontname','helvetica')
 

if range == [2 inf],
   titletext = ['normalised trend ' filename ', no filtering'];
else
   titletext = ['normalised trend ' filename ', [' num2str(range) '] samples per cycle'];
end
title(titletext,'fontsize',20,'fontname','helvetica')
yyaxis right
plot(1:length(clusters),clusters,'r','linewidth',1);
orient landscape  % optimize for printing

%--------------------------------------------------------
% PLOT AUTOCOVARIANCES
figure, clf reset
subplot(1,2,1), hold on
offset = 1.2;
for i = 1:NC,
   plot(autocorrmat(1:maxautocorr,i)-offset*i,'k','linewidth',2)
end

axis([0 maxautocorr -(NC+1)*offset 0]),grid
set(gca,'YTick',[-NC*offset:offset:-offset])
set(gca,'Yticklabel',tagnames(NC:-1:1))
set(gca,'fontsize',14,'fontname','helvetica','linewidth',1.5,'box','on')
xlabel('lag/sample interval','fontsize',14,'fontname','helvetica')
title(['autocovariance ' filename],'fontsize',14,...
   'fontname','helvetica')

%--------------------------------------------------------
% PLOT SPECTRA

% For plotting, fix the max value
% plotspecmat = 1.1*plotspecmat/max(max(plotspecmat));

% Axis for the plots. Top frequency is one half of 
% sampling frequency
xax1=([0:1/numfreq:(1-1/numfreq)])/2;% A half channel offset is not needed
% because it is not plotted as a bar spectrum. The spectrum starts at f = 0. 
index1 = 1:length(xax1);

xax=([1/(2*numfreq):1/(4*numfreq):(1-1/numfreq)])/2; % axis for splined spectra
index2 = 1:length(xax);

% replace spectra by splined curves for better visualisation
YY = spline(xax1(index1),plotspecmat(index1,:)',xax(index2));
plotspecmat = YY';

subplot(1,2,2)

offset = 1.2;
maxfreq = 1;  % axis is in fractions of sampling frequency

for i = 1:NC,
   %fix negative values(these are artefacts from splining)
   indexzero = find(plotspecmat(:,i)<0);
   plotspecmat(indexzero,i) = 0;
%    semilogx(xax,plotspecmat(:,i)-offset*i,'k','linewidth',2), hold on
   plot(xax,plotspecmat(:,i)-offset*i,'k','linewidth',2), hold on

end

tagnums = cellstr(num2str((1:NC)'));
axis([0 maxfreq -(NC+1)*offset 0]),grid on
set(gca,'YTick',[-NC*offset:offset:-offset])
set(gca,'Yticklabel',tagnums(NC:-1:1))
set(gca,'Xtick',[0;0.1;0.2;0.3;0.4;0.5])
set(gca,'fontsize',14,'fontname','helvetica','linewidth',1.5,'box','on')
xlabel('frequency/sampling frequency','fontsize',14,'fontname','helvetica')

title(['spectra ' filename], 'fontsize',14,'fontname','helvetica')
orient landscape  % optimize for printing

% figure(1)  % bring the time trend plot to the front

%=====================================================================================