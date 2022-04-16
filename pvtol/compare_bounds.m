clear;clf;
figure(1)
load("L1_ad_deccm\1e_4\CDC_RD_CCM_T_0.0001_lam_0.8_w_dist_1bound0.1_00_810_w_obs.mat");
% h1=plot(energyTraj);
ccm_dist = uTraj(1,:)-estDistTraj(1,:);
ccm = uTraj(1,:);
h1=plot(ccm);

hold on
load("L1_ad_deccm\1e_4\CDC_RD_CCM_T_0.0001_lam_0.8_with_dist_1bound0.1_00_810_w_obs.mat");
% h2=plot(energyTraj);
% h2=plot(uTraj(1,:));

load("L1_ad_deccm\1e_4\sim_L1_CCM_fixed_T_0.0001_lam_0.8_w_dist_1bound0.1_00_810_w_obs.mat");
% h3=plot(energyTraj);
% h3=plot(uTraj(1,:));


legend([h1,h2,h3],{'without','with', 'L1'})

% print('L1_ad_deccm/1e_4/compare_inputs_bound.pdf', '-painters', '-dpdf', '-r150');

figure(2)
load("L1_ad_deccm\1e_4\CDC_RD_CCM_T_0.0001_lam_0.8_w_dist_1bound0.1_00_810_w_obs.mat");
h4=plot(energyTraj);

hold on
load("L1_ad_deccm\1e_4\CDC_RD_CCM_T_0.0001_lam_0.8_with_dist_1bound0.1_00_810_w_obs.mat");
h5=plot(energyTraj);

load("L1_ad_deccm\1e_4\sim_L1_CCM_fixed_T_0.0001_lam_0.8_w_dist_1bound0.1_00_810_w_obs.mat");
h6=plot(energyTraj);

legend([h4,h5,h6],{'without','with', 'L1'})