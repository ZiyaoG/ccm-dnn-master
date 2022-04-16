load("L1_ad_deccm\1e_4\CDC_RD_CCM_T_0.0001_lam_0.8_w_dist_1bound0.1_00_810_w_obs.mat");
h1=plot(energyTraj);
hold on
load("L1_ad_deccm\1e_4\CDC_RD_CCM_T_0.0001_lam_0.8_with_dist_1bound0.1_00_810_w_obs.mat");
h2=plot(energyTraj);
load("L1_ad_deccm\1e_4\sim_L1_CCM_fixed_T_0.0001_lam_0.8_w_dist_1bound0.1_00_810_w_obs.mat");
h3=plot(energyTraj);

legend([h1,h2,h3],{'without','with', 'L1'})