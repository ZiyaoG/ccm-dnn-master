% Calculate Energy
% calculate disturbance wrt time
clear;
close all;

% Task1: 00_810
load('sim_ccm_lam_0.8_w_dist_1_00_810_w_obs.mat')
R_energy_ccm_no = energyTraj;
E1 = plot(times,R_energy_ccm_no,'r-','Linewidth',1);
hold on

load('sim_de_ccm_lam_0.8_w_dist_1_00_810_w_obs.mat')
R_energy_deccm_no = energyTraj;
E2 = plot(times,R_energy_deccm_no,'b-','Linewidth',1);
load('sim_ccm_lam_0.8_w_dist_1_with_poor_Adam__00_810_w_obs.mat')
R_energy_ccm_poor = energyTraj;
E3 = plot(times,R_energy_ccm_poor,'r-.','Linewidth',1);
load('sim_de_ccm_lam_0.8_w_dist_1_with_poor_Adam__00_810_w_obs.mat')
R_energy_deccm_poor = energyTraj;
E4 = plot(times,R_energy_deccm_poor,'b-.','Linewidth',2);
load('sim_de_ccm_lam_0.8_w_dist_1_with_perfect_Adam__00_810_w_obs.mat')
R_energy_ccm_perfect = energyTraj;
E5 = plot(times,R_energy_ccm_perfect,'r:','Linewidth',2.5);
load('sim_de_ccm_lam_0.8_w_dist_1_with_perfect_Adam__00_810_w_obs.mat')
R_energy_deccm_perfect = energyTraj;
E6 = plot(times,R_energy_deccm_perfect,'b:','Linewidth',2);

axis square
xlabel('Time (s)','interpreter','latex')
ylabel('Energy','interpreter','latex')
legend([E1,E2,E3,E4,E5,E6],{'CCM: no L', 'DE-CCM: no L', 'CCM: poor L', 'DE-CCM: poor L', 'CCM: good L', 'DE-CCM: good L'},'NumColumns',2,'Location','northwest','Orientation','horizontal');

xlim([0 10.5]);
% ylim([0 12]);
goodplot([6 5]);

% print('With partial learning_ccm fails.pdf', '-painters', '-dpdf', '-r150');
print('Plot Renergy .pdf', '-painters', '-dpdf', '-r150');
