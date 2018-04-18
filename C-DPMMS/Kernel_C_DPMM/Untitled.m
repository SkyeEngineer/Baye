 figure('Name','Fault Detection Results');

plot(1:length(rs),rs,'-r','Linewidth',3);
xlabel('Time(s)');ylabel('Results')
names = {'Normal'; 'Abnormal';};
set(gca,'ytick',[0,1],'yticklabel',names)
ylim([0-0.5 1+0.5]);
hold on; plot([0 0],[-0.5 1.5],'--b','Linewidth',1)
hold on; plot([901 901],[-0.5 1.5],'--k','Linewidth',1)
hold on; plot([3661 3661],[-0.5 1.5],'--b','Linewidth',1)
hold on; plot([4802 4802],[-0.5 1.5],'--k','Linewidth',1)
hold on; plot([11504 11504],[-0.5 1.5],'--b','Linewidth',1)
hold on; plot([11805 11805],[-0.5 1.5],'--k','Linewidth',1)
