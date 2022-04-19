% calculate disturbance wrt time
clear;
close all;

% Task1: 00_810
load('theSwedenWorkshop/data/_dtx_DECCM_ode1_constantintensity_03_03_00_88.mat');
dist_true = zeros(2,length(times));
for i=1:length(times)
    dist_true(:,i) = actual_dist_fcn(times(1,i),xTraj(:,i),[4,4],4);
end
hold on
h1 = plot(times,dist_true(1,:),'k--','Linewidth',1);
h2 = plot(times,dist_true(2,:),'k-.','Linewidth',1);
h3 = plot(times,estDistTraj(1,:),'r--','Linewidth',1);
h4 = plot(times,estDistTraj(2,:),'r-.','Linewidth',1);


axis square
xlabel('Time (s)','interpreter','latex')
ylabel('Disturbance','interpreter','latex')
legend([h1,h2,h3,h4],{'Actual 1', 'Actual 2', 'Estimated 1', 'Estimated 2'},'NumColumns',2,'Location','north','Orientation','vertical','interpreter','latex');

xlim([0 10.5]);
% ylim([-0.7 0.05]);
goodplot([5 5]);

% print('With partial learning_ccm fails.pdf', '-painters', '-dpdf', '-r150');
% print('Plot disturbance_test.pdf', '-painters', '-dpdf', '-r150');
%% some functions


% ---------------- True disturbance --------------------
function dist_force = actual_dist_fcn(t,x,center,radius)
% compute the disturbance force given x
% x is a n by m matrix, where each column represents a state vector value 
max_damping = 0.5;

% [dist_intensity,~] = dist_distribution(x(1,:),x(2,:),center,radius);
dist_intensity = 0.3;
dist_force_max = (x(4,:)^2+x(5,:)^2)*max_damping;
dist_force = [-1+0.3*sin(2*t); -1+0.3*cos(2*t)]*(dist_intensity.*dist_force_max); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% changed %%%%%%%%%%%%%%
% dist_force = -1./(1+exp(-5*phi))*0.1*(x(4)^2+x(5)^2);
end