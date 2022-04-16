% Calculate Energy
% calculate disturbance wrt time
clear;
close all;
figure(2)
% Task1: 00_810
load('L1_ad_deccm/1e_4/CDC_RD_CCM_T_0.0001_lam_0.8_w_dist_1bound0.1_00_810_w_obs.mat');
deccm = uTraj(1,:);
E1 = plot(times,deccm,'r-','Linewidth',1);
hold on

load('L1_ad_deccm/1e_4/sim_ad_ccm_T_0.0001_lam_0.8_w_dist_1_00_810_w_obs.mat');
ad = uTraj(1,:);
E2 = plot(times,ad,'b-','Linewidth',1);
load('L1_ad_deccm/1e_4/sim_L1_CCM_fixed_T_0.0001_lam_0.8_w_dist_1bound0.1_00_810_w_obs.mat');
L1 = uTraj(1,:);
E3 = plot(times,L1,'r-.','Linewidth',1);

axis square
xlabel('Time (s)','interpreter','latex')
ylabel('Control Input','interpreter','latex')
legend([E1,E2,E3],{'DECCM', 'ad', 'L1'},'NumColumns',3,'Location','northwest','Orientation','horizontal');

xlim([0 10.5]);
% ylim([0 12]);
goodplot([5 5]);

% print('With partial learning_ccm fails.pdf', '-painters', '-dpdf', '-r150');
% print('L1_ad_deccm/1e_4/Control Input.pdf', '-painters', '-dpdf', '-r150');
