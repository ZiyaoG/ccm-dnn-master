clear;clf;
close all;

figure(1);
hold on;
%% Task 1
load('adaptiveCCM_ode1_0.001_no_learning_00_1010_realbasefunc.mat');

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
h2 = plot(xTraj(1,:),xTraj(2,:),'b-','Linewidth',1.5);


h11 = scatter(xTraj(1,1),xTraj(2,1),'ro');
h12 = scatter(xTraj(1,end),xTraj(2,end),60,'rp');
h13 = scatter(xnomTraj(1,1),xnomTraj(2,1),'k*');

% %% Task 2
% load('adaptiveCCM_ode1_0.001_no_learning_00_1010_realbasefunc.mat');
% 
% h5 = plot(xTraj(1,:),xTraj(2,:),'b:','Linewidth',1.5);
% % h4 = plot(xnomTraj(1,:),xnomTraj(2,:),'k--','Linewidth',1.5);
% 
% scatter(xTraj(1,1),xTraj(2,1),'ro')
% scatter(xTraj(1,end),xTraj(2,end),60,'rp')
% scatter(xnomTraj(1,1),xnomTraj(2,1),'k*')


%% Task 3
load('sim_de_ccm_T_0.0001_lam_0.8_w_dist_1bound0.1_00_810_w_obs.mat');

h8 = plot(xTraj(1,:),xTraj(2,:),'r-','Linewidth',1.5);
% h7 = plot(xnomTraj(1,:),xnomTraj(2,:),'k--','Linewidth',1.5);
%% Task 2
load('robust_CCM_nonoise_no_learning_02ms_00_810.mat');

h4 = plot(xTraj(1,:),xTraj(2,:),'m-.','Linewidth',1.5);

%% Task 4
load('sim_ccm_lam_0.8_w_dist_1_00_810_w_obs.mat');

h5 = plot(xTraj(1,:),xTraj(2,:),'c:','Linewidth',1.5);

%% planned
load('adaptiveCCM_ode1_0.001_no_learning_00_1010_realbasefunc.mat');
h1 = plot(xnomTraj(1,:),xnomTraj(2,:),'k--','Linewidth',1.5);


scatter(xTraj(1,1),xTraj(2,1),'ro')
scatter(xTraj(1,end),xTraj(2,end),60,'rp')
scatter(xnomTraj(1,1),xnomTraj(2,1),'k*')
axis square
xlabel('$p_x$ (m)','interpreter','latex')
ylabel('$p_z$ (m)','interpreter','latex')
goodplot([5 5]);

%%
axis square
xlabel('$p_x$ (m)','interpreter','latex')
ylabel('$p_z$ (m)','interpreter','latex')
% legend([h1,h2,h3,h4,h5,h6,h7,h8,h9],{'T1:Planned', 'T1:DE-CCM', 'T1:CCM', 'T2:Planned', 'T2:DE-CCM', 'T2:CCM', 'T3:Planned', 'T3:DE-CCM', 'T3:CCM'},'NumColumns',3,'Location','North','Orientation','horizontal');
legend([h1,h8,h2,h4,h5],{'Planned', 'DE-CCM (ours)', 'Ad-CCM (true regressors)',  'RCCM', 'CCM'},'NumColumns',1,'Location','northwest','Orientation','horizontal','interpreter','latex');

xlim([-0 10.5]);
ylim([0 10.5]);
goodplot([5 5]);

print('no learning_adccm_deccm.pdf', '-painters', '-dpdf', '-r150');
%% functions
function [intensity,distance_to_center] = dist_distribution(X,Z,center,radius)
distance_to_center = sqrt((X-center(1)).^2 + (Z-center(2)).^2);


% --------------- using an inverse function ----------------
intensity = 1./(distance_to_center.^2+1);

% --------------- using an inverse function 2 ----------------
% intensity = 1./(sqrt(distance_to_center.^2)+1);
end