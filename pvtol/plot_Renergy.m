% Calculate Energy
% calculate disturbance wrt time
clear;
clf;
% close all;
figure(1)
% Task1: 00_810
file1 = load('theSwedenWorkshop/data/_dtx_DECCM_ode1_constantintensity_03_03_00_88.mat');
R_energy_ccm_no = file1.energyTraj;
times1 = file1.times;
E1 = plot(times1,R_energy_ccm_no,'r-','Linewidth',1);
hold on


file2 = load('theSwedenWorkshop/data/_dtx_CCM_ode1_constantintensity_03_03_00_88.mat');
R_energy_deccm_no = file2.energyTraj;
times2 = file2.times;
E2 = plot(times2,R_energy_deccm_no,'m:','Linewidth',1);

% file3 = load('simulation_results/sim_ccm_T_0.0001_lam_0.8_w_dist_1_with_poor_Adam_bound0.1_00_810_w_obs.mat');
file3 = load('theSwedenWorkshop/data/_dtx_adaptiveCCM_ode1_constantintensity_03_03_00_88.mat');
R_energy_ccm_poor = file3.energyTraj;
times3 = file3.times;
E3 = plot(times3,R_energy_ccm_poor,'b-.','Linewidth',1);

% file4 = load('simulation_results/sim_de_ccm_T_0.0001_lam_0.8_w_dist_1_with_poor_Adam_bound0.1_00_810_w_obs (1).mat');
file4 = load('theSwedenWorkshop/data/_dtx_adaptiveCCM_ode23_constantintensity_03_03_00_88_gain100_tVec.mat');
R_energy_deccm_poor = file4.energyTraj;
% times4 = file4.times;
E4 = plot(file4.tVec,R_energy_deccm_poor,'g-.','Linewidth',1);
% 
% file5 = load('simulation_results/sim_ccm_T_0.0001_lam_0.8_w_dist_1_with_perfect_Adam_bound0.1_00_810_w_obs.mat');
% R_energy_ccm_perfect = file5.energyTraj;
% times5 = file5.times;
% E5 = plot(times5,R_energy_ccm_perfect,'r:','Linewidth',2.5);
% 
% file6 = load('simulation_results/sim_de_ccm_T_0.0001_lam_0.8_w_dist_1_with_perfect_Adam_bound0.1_00_810_w_obs.mat');
% R_energy_deccm_perfect = file6.energyTraj;
% times6 = file6.times;
% E6 = plot(times6,R_energy_deccm_perfect,'b:','Linewidth',2);



axis square
xlabel('Time (s)','interpreter','latex')
ylabel('Energy','interpreter','latex')
% legend([E1,E2,E3,E4,E5,E6],{'CCM: no L', 'DE-CCM: no L', 'CCM: moderate L', 'DE-CCM: moderate L', 'CCM: good L', 'DE-CCM: good L'},'NumColumns',2,'Location','north','Orientation','horizontal');
legend([E1,E2,E3,E4],{'DECCM', 'CCM','AdCCM($\Gamma=10$)','AdCCM($\Gamma=100$)'},'NumColumns',2,'Location','north','Orientation','horizontal','interpreter','latex');

xlim([0 10]);
ylim([0 35]);
goodplot([5 5]);

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


% print('theSwedenWorkshop/Renergy_deccm_ad_ccm.pdf', '-painters', '-dpdf', '-r150');
