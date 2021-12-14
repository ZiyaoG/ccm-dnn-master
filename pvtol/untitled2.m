load('sim_de_ccm_2__T_0.0001_lam_0.8_w_dist_1_with_poor_Adam_bound0.1_00_810_w_obs.mat')
figure(1)
plot(times,xTraj(6,:))
% hold on
% plot(times,xTraj(5,:))
figure(2)
% plot(xTraj(1,:),xTraj(2,:))
plot(times,energyTraj)