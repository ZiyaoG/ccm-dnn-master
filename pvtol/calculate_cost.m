clear;
close all;

% Task1: 00_810
% Task2: 100_210
% Task3: 04_106
load('sim_ccm_lam_0.8_w_dist_1_with_poor_Adam__100_210_w_obs.mat')

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
%% results
y = [188.1646 185.6320 189.0798; 192.4614 206.9460 180.8110; 179.0480 202.9205 185.9694
; 178.6475 182.9514 173.2881
];
h = bar(y);
ylim([170 210])
name = {'No Learning';'Poor Learning';'Moderate Learning';'Good Learning'};
set(gca,'xticklabel',name)
ylabel('Cost')
set(h, {'DisplayName'}, {'Task 1','Task 2','Task 3'}')
% Legend will show names for each color
legend() 
goodplot([6 5])
print('Comparing Cost.pdf', '-painters', '-dpdf', '-r150');