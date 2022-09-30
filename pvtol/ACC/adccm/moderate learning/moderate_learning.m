clear;clf;
close all;

figure(1);
hold on;
%% Task 1
load('DECCM_moderate_learning_02ms_00_810.mat');

xx = -1.5:0.05:11;
zz = 0:0.05:11;
[X,Z] = meshgrid(xx,zz);
Dist_intensity= dist_distribution(X,Z,dist_config.center,dist_config.radius);
Dist_distribution.X = X;
Dist_distribution.Z = Z;
Dist_distribution.intensity = Dist_intensity;

obs = [3.1 6 0.6;           
    6.9 6 0.6;
    5 3 0.5];
gray_color = [1 1 1]*80/255;
trajGen_config.include_obs = 1;
trajGen_config.obs = obs;


if sim_config.include_dist == 1
    visualize_dist_area(Dist_distribution);
end
if trajGen_config.include_obs == 1
    visualize_obs(trajGen_config.obs,gray_color);
end
h2 = plot(xTraj(1,:),xTraj(2,:),'r-','Linewidth',1.5);


h11 = scatter(xTraj(1,1),xTraj(2,1),'ro');
h12 = scatter(xTraj(1,end),xTraj(2,end),60,'rp');
h13 = scatter(xnomTraj(1,1),xnomTraj(2,1),'k*');


%% Task 3
load('adaptiveCCM_ode1_0.0002_moderate_learning_00_810_polybasefunc_gain01.mat');

h8 = plot(xTraj(1,:),xTraj(2,:),'b-','Linewidth',1.5);
% h7 = plot(xnomTraj(1,:),xnomTraj(2,:),'k--','Linewidth',1.5);
%% Task 2
load('robust_CCM_moderate_learning_02ms_00_810.mat');

h4 = plot(xTraj(1,:),xTraj(2,:),'m-.','Linewidth',1.5);
%% planned
load('DECCM_moderate_learning_02ms_00_810.mat');
h1 = plot(xnomTraj(1,:),xnomTraj(2,:),'k--','Linewidth',1.5);

scatter(xTraj(1,1),xTraj(2,1),'ro')
scatter(xTraj(1,end),xTraj(2,end),60,'rp')
scatter(xnomTraj(1,1),xnomTraj(2,1),'k*')
axis square
xlabel('$p_x$ (m)','interpreter','latex')
ylabel('$p_z$ (m)','interpreter','latex')
goodplot([5 5]);
%% plot the trajactory with nn

% figure(2);
% hold on;
% 
% load('sim_deccm_lam_0.8_w_dist_1_60_7.510_w_obs.mat');
% 
% xx = -2.5:0.05:10.5;
% zz = 0:0.05:13;
% [X,Z] = meshgrid(xx,zz);
% Dist_intensity= dist_distribution(X,Z,dist_config.center,dist_config.radius);
% Dist_distribution.X = X;
% Dist_distribution.Z = Z;
% Dist_distribution.intensity = Dist_intensity;
% 
% obs = [3.4 6 0.6;           
%         6.6 6 0.6];
% gray_color = [1 1 1]*80/255;
% trajGen_config.include_obs = 1;
% trajGen_config.obs = obs;
% 
% 
% if sim_config.include_dist == 1
%     visualize_dist_area(Dist_distribution);
% end
% if trajGen_config.include_obs == 1
%     visualize_obs(trajGen_config.obs,gray_color);
% end
% h2 = plot(xTraj(1,:),xTraj(2,:),'r-','Linewidth',1.5);
% h1 = plot(xnomTraj(1,:),xnomTraj(2,:),'k:','Linewidth',1.5);
% scatter(xnomTraj(1,1),xnomTraj(2,1),'ko')
% scatter(xnomTraj(1,end),xnomTraj(2,end),'kp')
% 
% load('sim_ccm_lam_0.8_w_dist_1_withnn__60_7.510_w_obs.mat');
% h3 = plot(xTraj(1,:),xTraj(2,:),'b-','Linewidth',1.5);
% 
% 
% load('sim_deccm_lam_0.8_w_dist_1_withnn__00_1010_w_obs.mat');
% h5 = plot(xTraj(1,:),xTraj(2,:),'r-.','Linewidth',1.5);
% h4 = plot(xnomTraj(1,:),xnomTraj(2,:),'k--','Linewidth',1.5);
% scatter(xnomTraj(1,1),xnomTraj(2,1),'ko')
% scatter(xnomTraj(1,end),xnomTraj(2,end),'kp')
% 
% load('sim_ccm_lam_0.8_w_dist_1_withnn__00_1010_w_obs.mat');
% h6 = plot(xTraj(1,:),xTraj(2,:),'b-.','Linewidth',1.5);


%%
axis square
xlabel('$p_x$ (m)','interpreter','latex')
ylabel('$p_z$ (m)','interpreter','latex')
% legend([h1,h2,h3,h4,h5,h6,h7,h8,h9],{'T1:Planned', 'T1:DE-CCM', 'T1:CCM', 'T2:Planned', 'T2:DE-CCM', 'T2:CCM', 'T3:Planned', 'T3:DE-CCM', 'T3:CCM'},'NumColumns',3,'Location','North','Orientation','horizontal');
legend([h1, h2,h8,h4],{'Planned', 'DE-CCM (ours)', 'Ad-CCM (poly. regressors)',  'RCCM'},'NumColumns',1,'Location','northwest','Orientation','horizontal','interpreter','latex');

xlim([-0 10.5]);
ylim([0 10.5]);
goodplot([5 5]);

print('moderate learning_adccm_deccm.pdf', '-painters', '-dpdf', '-r150');
%% functions
function [intensity,distance_to_center] = dist_distribution(X,Z,center,radius)
distance_to_center = sqrt((X-center(1)).^2 + (Z-center(2)).^2);


% --------------- using an inverse function ----------------
intensity = 1./(distance_to_center.^2+1);

% --------------- using an inverse function 2 ----------------
% intensity = 1./(sqrt(distance_to_center.^2)+1);
end