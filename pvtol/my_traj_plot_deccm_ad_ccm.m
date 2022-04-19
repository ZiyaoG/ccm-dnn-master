clear;clf;
close all;

figure(1);
hold on;
%% Task 1

load('theSwedenWorkshop/data/_dtx_CCM_ode1_constantintensity_03_03_00_88.mat');
xx = -1.5:0.05:13;
zz = 0:0.05:13;
[X,Z] = meshgrid(xx,zz);
Dist_intensity= dist_distribution(X,Z,dist_config.center,dist_config.radius);
Dist_distribution.X = X;
Dist_distribution.Z = Z;
Dist_distribution.intensity = Dist_intensity;

obs = [2.1 5 0.6;           
    5.9 5 0.6;
    4 2 0.6];
gray_color = [1 1 1]*80/255;
trajGen_config.include_obs = 1;
trajGen_config.obs = obs;


% if sim_config.include_dist == 1
%     visualize_dist_area(Dist_distribution);
% end
if trajGen_config.include_obs == 1
    visualize_obs(trajGen_config.obs,gray_color);
end
h1 = plot(xnomTraj(1,:),xnomTraj(2,:),'k--','Linewidth',3);

h3 = plot(xTraj(1,:),xTraj(2,:),'Color', 'm','LineStyle',':','Linewidth',1.5);

% h11 = scatter(xTraj(1,1),xTraj(2,1),'ro');
% h12 = scatter(xTraj(1,end),xTraj(2,end),60,'rp');
% h13 = scatter(xnomTraj(1,1),xnomTraj(2,1),'k*');
%% 
load('theSwedenWorkshop/data/_dtx_adaptiveCCM_ode23_constantintensity_03_03_00_88_gain100.mat');
h5 = plot(xTraj(1,:),xTraj(2,:),'Color', [0 1 0],'LineStyle','-.','Linewidth',2);
%% 
load('theSwedenWorkshop/data/_dtx_adaptiveCCM_ode1_constantintensity_03_03_00_88.mat');
h4 = plot(xTraj(1,:),xTraj(2,:),'Color', [0 0 1],'LineStyle','-.','Linewidth',1);
% [0.4660 0.6740 0.1880][0 0.4470 0.7410]

%% 
load('theSwedenWorkshop/data/_dtx_DECCM_ode1_constantintensity_03_03_00_88.mat');
h2 = plot(xTraj(1,:),xTraj(2,:),'Color', 'r','LineStyle','-','Linewidth',1);
%%
axis square
xlabel('$p_x$ (m)','interpreter','latex')
ylabel('$p_z$ (m)','interpreter','latex')
legend([h1,h2,h3,h4,h5],{'Planned', 'DE-CCM','CCM', 'Ad-CCM($\Gamma=10$)','Ad-CCM($\Gamma=100$)'},'NumColumns',3,'Location','Northwest','Orientation','horizontal','interpreter','latex');

xlim([0 8.5]);
ylim([0 8.5]);
goodplot([5 5]);

print('theSwedenWorkshop/Traj_deccm_ad_ccm.pdf', '-painters', '-dpdf', '-r150');
%% functions
function [intensity,distance_to_center] = dist_distribution(X,Z,center,radius)
distance_to_center = sqrt((X-center(1)).^2 + (Z-center(2)).^2);


% --------------- using an inverse function ----------------
intensity = 1./(distance_to_center.^2+1);

% --------------- using an inverse function 2 ----------------
% intensity = 1./(sqrt(distance_to_center.^2)+1);
end