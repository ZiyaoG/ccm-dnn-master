% Calculate Energy
% calculate disturbance wrt time
clear;
close all;
figure(2)
% Task1: 00_810

load('robust_CCM_nonoise_no_learning_02ms_delay50steps_00_810.mat');
adccm1 = uTraj(1,:);
E3 = plot(times,adccm1,'b-.','Linewidth',0.7);
hold on
load('DECCM_no_learning_02ms_delay50steps_delayinestimator_00_810.mat');
adccm23 = uTraj(1,:);
E4 = plot(times,adccm23,'-r','Linewidth',1);


% axis square
xlabel('Time (s)','interpreter','latex')
ylabel('$u_1$ ($N$)','interpreter','latex')
legend([E3,E4],{'RCCM', 'DE-CCM (ours)'},'NumColumns',2,'Location','north','Orientation','horizontal','interpreter','latex');

xlim([0 10.5]);
% ylim([0 12]);
goodplot([5 5]);

ax = gca;
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
set(ax,'ticklabelinterpreter','none') %or 'tex' but not 'latex'

% print('deccm_Rccm_inputs_w_delays.pdf', '-painters', '-dpdf', '-r150');
% print('theSwedenWorkshop/Input1_deccm_ad_ccm.pdf', '-painters', '-dpdf', '-r150');