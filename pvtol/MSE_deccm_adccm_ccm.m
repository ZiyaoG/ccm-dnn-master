clear;clf;
close all;
dim = 2;
figure(1);
%% Task 1
figure(1)
load('theSwedenWorkshop/data/_dtx_DECCM_ode1_constantintensity_03_03_00_88.mat');
% subplot(3,1,1)
hold on;


delta_x = xnomTraj-xTraj;
x_xnom = vecnorm(delta_x);
R = sqrt(controller.w_upper/controller.w_lower);
h1 = fplot(@(t) R*norm([0 0 0] - [5/180*pi 5/180*pi 0])*exp(-controller.lambda*t),'Color', 'k','LineStyle','--','Linewidth',1);

h2 = plot(times(1,1:end-1),x_xnom(1,1:end-1),'Color', 'r','LineStyle','-','Linewidth',1);
xlim([0 13]);
MSE_deccm1 = norm(delta_x,"fro")^2/100001;
% MSE_deccm2 = mean(x_xnom.^2);
Int_deccm = sum(x_xnom.^2*1e-4)/10;
%% 
load('theSwedenWorkshop/data/_dtx_adaptiveCCM_ode23_constantintensity_03_03_00_88_gain100_tVec.mat');
% subplot(3,1,1)
load('nomTraj_x0xF_config_1_no_learning.mat')
x_nom_fcn = soln.interp.state;
u_nom_fcn = soln.interp.control;
times = tVec';
simuLen = length(times);
xnomTraj = zeros(6,simuLen);
unomTraj = zeros(2,simuLen);
for t =1:simuLen-1
    xnomTraj(:,t) = x_nom_fcn(times(t));
    unomTraj(:,t) = u_nom_fcn(times(t));
end
delta_x = xnomTraj-xTraj;
x_xnom = vecnorm(delta_x);
% times = times';
h7 = plot(times(1,1:end-1),x_xnom(1,1:end-1),'Color', 'b','LineStyle','-.','Linewidth',1);
MSE_ad100 = norm(delta_x,"fro")^2/4044;
% MSE_adccm2 = mean(x_xnom.^2);
% adx = x_xnom;
dt = times(1,2:end) - times(1,1:end-1);
Int_adccm1 = sum(dt.*(x_xnom(2:end).^2))/10;
Int_adccm2 = sum(dt.*(x_xnom(1:end-1).^2))/10;
Int_adccm = (Int_adccm1+Int_adccm2)/2;
%% 
load('theSwedenWorkshop/data/_dtx_CCM_ode1_constantintensity_03_03_00_88.mat');
% subplot(3,1,1)
delta_x = xnomTraj-xTraj;
x_xnom = vecnorm(delta_x);
h10 = plot(times(1,1:end-1),x_xnom(1,1:end-1),'Color', 'm','LineStyle',':','Linewidth',1);
MSE_ccm = norm(delta_x,"fro")^2/100001;
% MSE_ccm2 = mean(x_xnom.^2);
% ccmx = x_xnom;
Int_ccm = sum(x_xnom.^2*1e-4)/10;
%%
% subplot(3,1,1)
xlabel('Time (s)','interpreter','latex')
ylabel('$\|x-x^*\|$','interpreter','latex')
legend([h2,h7,h10],{'DE-CCM','Ad-CCM', 'CCM'},'NumColumns',3,'Location','North','Orientation','horizontal','interpreter','latex');
xlim([0 13]);
% ylim([0 0.23]);
goodplot([6.5 6.5]);

% print('figures/x_xnom_deccm_ad_ccm.pdf', '-painters', '-dpdf', '-r150');

%% calculate MSE

% figure(2)
% mean(ccmx)
% mean(adx)
