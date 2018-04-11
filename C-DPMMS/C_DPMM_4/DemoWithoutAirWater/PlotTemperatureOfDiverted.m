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


Training=[Normal_120air_01water' Normal_150air_05water' ];
figure(1),
subplot(2,1,1)
l1=plot(Diverted_120air_01water(:,3),'Linewidth',3);
hold on,l2=plot(Diverted_150air_05water(:,3),'Linewidth',3);
hold on ,l3=plot(Training(3,:),'Linewidth',3);

hold on; %plot(601,901,541,1141,301)
a=601;b=a+901;c=b+541;d=c+1141;e=d+301;
plot([a a],ylim,'LineStyle','--');text(a-350,19,'day1')
hold on
line([b b],ylim,'LineStyle','--');text(b-350,19,'day2')
hold on
line([c c],ylim,'LineStyle','--');text(c-350,19,'day1')
hold on
line([d d],ylim,'LineStyle','--');text(d-350,19,'day2')
hold on
line([e e],ylim,'LineStyle','--');text(e-300,19,'day2')
title('Air Temperature' )
names = {'Normal(120air,01water)-(A)'; 'Normal(150air,05water)-(B)'};
set(gca,'xtick',[length(Normal_120air_01water),length(Training)],'xticklabel',names)
legend([l1 l2 l3],{'Diverted(120air,01water)-(A)','Diverted(150air,05water)-(B)','Training'},'Location','best')
ylabel('Temperature')


subplot(2,1,2)
l4=plot(Diverted_120air_01water(:,7),'Linewidth',3);
hold on,l5=plot(Diverted_150air_05water(:,7),'Linewidth',3);
hold on ,l6=plot(Training(7,:),'Linewidth',3);
line([a a],ylim,'LineStyle','--');text(a-350,26,'day1')
hold on
line([b b],ylim,'LineStyle','--');text(b-350,26,'day2')
hold on
line([c c],ylim,'LineStyle','--');text(c-350,26,'day1')
hold on
line([d d],ylim,'LineStyle','--');text(d-350,26,'day2')
hold on
line([e e],ylim,'LineStyle','--');text(e-300,26,'day2')
title('Water Temperature')
names = {'Normal(120air,01water)-(A)'; 'Normal(150air,05water)-(B)'};
set(gca,'xtick',[length(Normal_120air_01water),length(Training)],'xticklabel',names)
legend([l4 l5 l6],{'Diverted(120air,01water)-(A)','Diverted{150air,05water)-(B)','Training'},'Location','best')
ylabel('Temperature')
% % print('TemperaturePlot','-depsc')





