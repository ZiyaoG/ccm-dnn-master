% Calculate Energy
% calculate disturbance wrt time
clear;
clf;
% close all;
figure(1)
% Task1: 00_810
file1 = load('ACC/RCCM_no_learning_1ms_0.03s_delay_00_810.mat');
% file1 = load('ACC/DECCM_no_learning_1ms_0.03s_delay_00_810.mat');

R_energy_ccm_no = file1.energyTraj;
times1 = file1.times;
E1 = plot(times1,R_energy_ccm_no,'b-.','Linewidth',1);
hold on


file2 = load('ACC/RCCM_perfect_learning_1ms_0.03s_delay_00_810.mat');
% file2 = load('ACC/DECCM_perfect_learning_1ms_0.03s_delay_00_810.mat');

R_energy_deccm_no = file2.energyTraj;
times2 = file2.times;
E2 = plot(times2,R_energy_deccm_no,'k:','Linewidth',1.5);

file3 = load('ACC/DECCM_no_learning_1ms_0.03s_delay_00_810.mat');
R_energy_ccm_poor = file3.energyTraj;
times3 = file3.times;
E3 = plot(times3,R_energy_ccm_poor,'m-','Linewidth',1.5);

file4 = load('ACC/DECCM_perfect_learning_1ms_0.03s_delay_00_810.mat');
R_energy_deccm_poor = file4.energyTraj;
times4 = file4.times;
E4 = plot(times4,R_energy_deccm_poor,'r-','Linewidth',2);0


% axis square
xlabel('Time (s)','interpreter','latex')
ylabel('Energy','interpreter','latex')
% legend([E1,E2,E3,E4,E5,E6],{'CCM: no L', 'DE-CCM: no L', 'CCM: moderate L', 'DE-CCM: moderate L', 'CCM: good L', 'DE-CCM: good L'},'NumColumns',2,'Location','north','Orientation','horizontal');
% legend([E1,E2,E3,E4],{'RCCM - no learning', 'RCCM - good learning','DE-CCM - no learning','DE-CCM - good learning'},'NumColumns',2,'Location','north','Orientation','horizontal','interpreter','latex');
% legend([E1,E2],{'DECCM - no learning', 'DECCM - good learning'},'NumColumns',2,'Location','north','Orientation','horizontal','interpreter','latex');

xlim([6 6.5]);
% ylim([0 35]);
goodplot([6.3 3]);

%%
% axes('position',[.3 .45 .5 .25])
% box on % put box around new pair of axes
% indexOfInterest = (times2 < 9.5) & (times2 > 4.5); % range of t near perturbation
% plot(times2(indexOfInterest),R_energy_deccm_no(indexOfInterest),'b-','Linewidth',1) % plot on new axes
% hold on
% indexOfInterest = (times4 < 9.5) & (times4 > 4.5);
% plot(times4(indexOfInterest),R_energy_deccm_poor(indexOfInterest),'b-.','Linewidth',2)
% indexOfInterest = (times5 < 9.5) & (times5 > 4.5);
% plot(times5(indexOfInterest),R_energy_ccm_perfect(indexOfInterest),'r:','Linewidth',2.5)
% indexOfInterest = (times6 < 9.5) & (times6 > 4.5);
% plot(times6(indexOfInterest),R_energy_deccm_perfect(indexOfInterest),'b:','Linewidth',2)
% ax = gca;
% ax.FontSize = 13;
% set(gca,'XTick',[5,7,9])
% set(gca,'XTickLabel',{'5','7','9'})
% set(gca,'YTick',[0,0.01,0.02])
% set(gca,'YTickLabel',{'0','0.01','0.02'})
% axis tight


print('ACC/Renergy_deccm_Rccm_30msdelay_zoomin.pdf', '-painters', '-dpdf', '-r150');
