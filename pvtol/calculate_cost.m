clear;
close all;

% Task1: 00_810 179.0226
% Task2: 100_210 190.7433
% Task3: 04_106 181.6375
load('simulation_results_CDC/safe_exploration/CDC_RD_CCM_T_0.0001_lam_0.8_w_dist_1_with_XX_Adam_bound0.1_100_210_w_obs.mat')

dt = sim_config.step_size;
u_size = size(uTraj);
cost = 0.0;
for i = 1:u_size(2)
    u_t = uTraj(:,i);
    cost_t = u_t(1)^2+u_t(2)^2;
    cost = cost + cost_t*dt;
end
cost = cost + 5*u_size(2)*dt;
disp(cost)
% return
%% results
% y = [188.1646 185.6320 189.0798; 192.4614 206.9460 180.8110; 179.0226 190.7433 181.6375
% ; 178.6475 182.9514 173.2881];
y = [188.1646 185.6320 189.0798;179.0226 190.7433 181.6375
; 178.6475 182.9514 173.2881];
h = bar(y);
ylim([170 195])
name = {'No Learning';'Moderate Learning';'Good Learning'};
set(gca,'xticklabel',name)
ylabel('Cost')
set(h, {'DisplayName'}, {'Task 1','Task 2','Task 3'}')
% Legend will show names for each color
legend() 
goodplot([6 5])
% print('CDC_figure/Comparing Cost_safe exploration.pdf', '-painters', '-dpdf', '-r150');