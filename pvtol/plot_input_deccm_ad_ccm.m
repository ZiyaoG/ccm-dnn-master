% Calculate Energy
% calculate disturbance wrt time
clear;
close all;
figure(2)
% Task1: 00_810

load('theSwedenWorkshop/data/_dtx_adaptiveCCM_ode1_constantintensity_03_03_00_88.mat');
adccm1 = uTraj(1,:);
E3 = plot(times,adccm1,'b-.','Linewidth',0.7);
hold on
load('theSwedenWorkshop/data/_dtx_adaptiveCCM_ode23_constantintensity_03_03_00_88_gain100_tVec.mat');
adccm23 = uTraj(1,:);
E4 = plot(tVec,adccm23,':g','Linewidth',0.8);


load('theSwedenWorkshop/data/_dtx_DECCM_ode1_constantintensity_03_03_00_88.mat');
deccm = uTraj(1,:);
E1 = plot(times,deccm,'r-','Linewidth',1.5);
% hold on

load('theSwedenWorkshop/data/_dtx_CCM_ode1_constantintensity_03_03_00_88.mat');
% load('theSwedenWorkshop/data/robust_CCM_ode1_constantintensity_03_03_00_88.mat');
ccm = uTraj(1,:);
E2 = plot(times,ccm,'m--','Linewidth',1.5);


axis square
xlabel('Time (s)','interpreter','latex')
ylabel('$u_1$ ($N$)','interpreter','latex')
legend([E1,E2,E3,E4],{'DE-CCM', 'CCM', 'Ad-CCM($\Gamma=10$)','Ad-CCM($\Gamma=100$)'},'NumColumns',2,'Location','north','Orientation','horizontal','interpreter','latex');

xlim([0 10.5]);
% ylim([0 12]);
goodplot([5 5]);

ax = gca;
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
set(ax,'ticklabelinterpreter','none') %or 'tex' but not 'latex'

% print('With partial learning_ccm fails.pdf', '-painters', '-dpdf', '-r150');
print('theSwedenWorkshop/Input1_deccm_ad_ccm.pdf', '-painters', '-dpdf', '-r150');
