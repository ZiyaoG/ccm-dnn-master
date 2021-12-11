% Calculate Energy
% calculate disturbance wrt time
clear;
close all;

% Task1: 00_810
load('simulation_results/sim_ccm_lam_0.8_w_dist_1_00_810_w_obs.mat')
Input_ccm_no = uTraj(1,:);
E1 = plot(times,Input_ccm_no,'r-','Linewidth',1);
hold on

load('simulation_results/sim_de_ccm_T_0.0001_lam_0.8_w_dist_1bound0.1_00_810_w_obs.mat')
Input_deccm_no = uTraj(1,:);
E2 = plot(times,Input_deccm_no,'b-','Linewidth',1);
load('simulation_results/sim_ccm_lam_0.8_w_dist_1_with_poor_Adam__00_810_w_obs.mat')
Input_ccm_poor = uTraj(1,:);
E3 = plot(times,Input_ccm_poor,'r-.','Linewidth',1);
load('simulation_results/sim_de_ccm_T_0.0001_lam_0.8_w_dist_1_with_poor_Adam_bound0.1_00_810_w_obs.mat')
Input_deccm_poor = uTraj(1,:);
E4 = plot(times,Input_deccm_poor,'b-.','Linewidth',1);
load('simulation_results/sim_ccm_lam_0.8_w_dist_1_with_perfect_Adam__00_810_w_obs.mat')
Input_ccm_perfect = uTraj(1,:);
E5 = plot(times,Input_ccm_perfect,'r:','Linewidth',1);
load('simulation_results/sim_de_ccm_lam_0.8_w_dist_1_with_perfect_Adam_bound0.1_00_810_w_obs.mat')
Input_deccm_perfect = uTraj(1,:);
E6 = plot(times,Input_deccm_perfect,'b:','Linewidth',1);

axis square
xlabel('Time (s)','interpreter','latex')
ylabel('Energy','interpreter','latex')
legend([E1,E2,E3,E4,E5,E6],{'CCM: no L', 'DE-CCM: no L', 'CCM: poor L', 'DE-CCM: poor L', 'CCM: good L', 'DE-CCM: good L'},'NumColumns',2,'Location','northwest','Orientation','horizontal');

xlim([0 10.5]);
% ylim([0 12]);
goodplot([6 5]);

% print('With partial learning_ccm fails.pdf', '-painters', '-dpdf', '-r150');
print('Plot input_bound0.1_T_0.0001.pdf', '-painters', '-dpdf', '-r150');
