n = plant.n; nu = plant.nu;
% ccm_sim_data = 'ccm_sim_w_dist_1.mat';
% rccm_sim_data = 'ccm_sim_w_dist_1.mat';
%-------------------------------------------------
% nominal trajectory
simuLen = length(times);
xnomTraj = zeros(n,simuLen);
unomTraj = zeros(nu,simuLen);
for t =1:simuLen
    xnomTraj(:,t) = x_nom_fcn(times(t));
    unomTraj(:,t) = u_nom_fcn(times(t));
end

addText = 'CCM';
close all;
figure(1); clf;
subplot(2,2,1); hold on;
plot(times,xnomTraj(1,:),'b--',times,xnomTraj(2,:),'r--');
plot(times,xTraj(1,:),'b-',times,xTraj(2,:),'r-');
xlabel('Time (s)')
ylabel('x & z (m)')
legend('x_{nom}', 'z_{nom}', ['x: ' addText],['z: ' addText]);

subplot(2,2,2); hold on;
plot(times,xnomTraj(3,:)*180/pi,'--');
plot(times,xTraj(3,:)*180/pi);
xlabel('Time (s)')
ylabel('\phi (deg)')
legend('Nominal',addText);

subplot(2,2,3);hold on;
plot(times,unomTraj(1,:),'b--',times,unomTraj(2,:),'r--');
plot(times,uTraj(1,:),'b-',times,uTraj(2,:),'r-');
xlabel('Time (s)');
ylabel('u (N)')
legend('u_{nom,1}', 'u_{nom,2}',['u_1: ' addText],['u_2: ' addText]);
% plot(times,distTraj);
% xlabel('Time (s)');
% ylabel('$\|x-x^\star\|/\|x_0-x^\star_0\|$','interpreter','latex');

subplot(2,2,4)
plot(times, energyTraj);
xlabel('Time (s)');
ylabel('Riemann energy')
legend(addText)

figure(2);
clf;hold on;
if sim_config.include_dist == 1
    plot(times,distTraj(1,:),'k',times,estDistTraj(1,:),'r--'); 
    plot(times,distTraj(2,:),'k--',times,estDistTraj(2,:),'r-.');
    ylabel('Disturbance')
    xlabel('Time (s)');
    legend('Actual: $\tilde{d}_1$','Estimated: $\hat{\sigma}_1$','Actual: $\tilde{d}_2$','Estimated: $\hat{\sigma}_2$','Interpreter','Latex');    
end

figure(3);clf;
hold on;
if sim_config.include_dist == 1
    visualize_dist_area(Dist_distribution);
end
if trajGen_config.include_obs == 1
    visualize_obs(trajGen_config.obs,gray_color);
end
h2 = plot(xTraj(1,:),xTraj(2,:),'r-','Linewidth',1.5);
h1 = plot(xnomTraj(1,:),xnomTraj(2,:),'k--','Linewidth',1.5);

scatter(xnomTraj(1,1),xnomTraj(2,1),'ko')
scatter(xnomTraj(1,end),xnomTraj(2,end),'k*')
scatter(xTraj(1,end),xTraj(2,end),'r*')
axis square
xlabel('p_x (m)')
ylabel('p_z (m)')
legend([h1,h2],{'Planned','RD-CCM'});
xlim([0 10]);
ylim([0 10]);
w_max = 1;
% if sim_config.save_sim_rst == 1
%     if controller.distEstScheme == 2    
%         file_name = 'CDC_RD_CCM';  
%     elseif controller.distEstScheme == 0
%         file_name = 'CDC_CCM';
%     end
%     file_name = [file_name '_T_' num2str(sim_config.step_size)];
%     file_name = [file_name '_lam_' num2str(controller.lambda,2)];
%     if sim_config.include_dist == 1
%         file_name = [file_name '_w_dist_' num2str(w_max)];
%     end
%     
%     if use_distModel_in_planning_control == 1
%         file_name = [file_name '_with_XX_Adam_'];
%     end
%     file_name = [file_name 'bound' num2str(distEst_config.est_errBnd)];
%     file_name = [file_name '_' num2str(x0(1)) num2str(x0(2)) '_' num2str(xF(1)) num2str(xF(2))];
%     if sim_config.include_obs == 1   
%         file_name  = [file_name '_w_obs.mat'];
%     else
%         file_name  = [file_name '.mat'];
%     end
%     save(file_name,'times','xTraj','uTraj','xnomTraj','unomTraj','energyTraj','dist_config','sim_config','plant','controller','estDistTraj');
% end