% Calculate Energy
% calculate disturbance wrt time
clear;
close all;
figure(2)
% Task1: 00_810
load('simulation_results/sim_ccm_T_0.0001_lam_0.8_w_dist_1bound0.1_00_810_w_obs.mat');
Input_ccm_no = uTraj(1,:);
E1 = plot(times,Input_ccm_no,'r-','Linewidth',1);
hold on

load('simulation_results/sim_de_ccm_T_0.0001_lam_0.8_w_dist_1bound0.1_00_810_w_obs.mat');
Input_deccm_no = uTraj(1,:);
E2 = plot(times,Input_deccm_no,'b-','Linewidth',1);
load('simulation_results_CDC/safe_exploration/CDC_CCM_T_0.0001_lam_0.8_w_dist_1_with_XX_Adam_bound0.1_00_810_w_obs.mat');
Input_ccm_poor = uTraj(1,:);
E3 = plot(times,Input_ccm_poor,'r-.','Linewidth',1);
load('simulation_results_CDC/safe_exploration/CDC_RD_CCM_T_0.0001_lam_0.8_w_dist_1_with_XX_Adam_bound0.1_00_810_w_obs.mat');
Input_deccm_poor = uTraj(1,:);
E4 = plot(times,Input_deccm_poor,'b-.','Linewidth',2);

load('simulation_results/sim_ccm_T_0.0001_lam_0.8_w_dist_1_with_perfect_Adam_bound0.1_00_810_w_obs.mat');
Input_ccm_good = uTraj(1,:);
E5 = plot(times,Input_ccm_good,'r:','Linewidth',2.5);
load('simulation_results/sim_de_ccm_T_0.0001_lam_0.8_w_dist_1_with_perfect_Adam_bound0.1_00_810_w_obs.mat');
Input_deccm_good = uTraj(1,:);
E6 = plot(times,Input_deccm_good,'b:','Linewidth',2);


axis square
xlabel('Time (s)','interpreter','latex')
ylabel('Control Input','interpreter','latex')
legend([E1,E2,E3,E4,E5,E6],{'CCM: no L', 'DE-CCM: no L', 'CCM: moderate L', 'DE-CCM: moderate L', 'CCM: good L', 'DE-CCM: good L'},'NumColumns',2,'Location','north','Orientation','horizontal');
% legend([E1,E2,E3,E4],{'DECCM min u', 'L1','DECCM min u-u^*','u planned'},'NumColumns',4,'Location','northwest','Orientation','horizontal');

xlim([0 10.5]);
ylim([0 5.3]);
goodplot([5 5]);

% print('With partial learning_ccm fails.pdf', '-painters', '-dpdf', '-r150');
print('CDC_figure/Control Input.pdf', '-painters', '-dpdf', '-r150');
