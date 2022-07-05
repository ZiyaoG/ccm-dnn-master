n = plant.n; nu = plant.nu;
% nominal trajectory
simuLen = length(times);
xnomTraj = zeros(n,simuLen);
unomTraj = zeros(nu,simuLen);
for t =1:simuLen
    xnomTraj(:,t) = x_nom_fcn(times(t));
    unomTraj(:,t) = u_nom_fcn(times(t));
end
times = 0:0.05:duration;
addText = 'CCM';
close all;
figure(1); clf;
subplot(2,2,1); hold on;
plot(times,xnomTraj(1,:),'b--',times,xnomTraj(2,:),'r--');
plot(times,xTraj(1,:),'b-',times,xTraj(2,:),'r-');
% plot(tVec,xTraj(1,:),'b-',tVec,xTraj(2,:),'r-');
xlabel('Time (s)')
ylabel('x & z (m)')
legend('x_{nom}', 'z_{nom}',['x: ' addText],['z: ' addText]);

subplot(2,2,2); hold on;
plot(times,xnomTraj(3,:)*180/pi,'--');
plot(times,xTraj(3,:)*180/pi);
% plot(tVec,xTraj(3,:)*180/pi);
xlabel('Time (s)')
ylabel('\phi (deg)')
legend('Nominal',addText);

subplot(2,2,3);hold on;
plot(times,unomTraj(1,:),'b--',times,unomTraj(2,:),'r--');
plot(times,uTraj(1,:),'b-',times,uTraj(2,:),'r-');
% plot(tVec,uTraj(1,:),'b-',tVec,uTraj(2,:),'r-');
xlabel('Time (s)');
ylabel('u (N)')
legend('u_{nom,1}', 'u_{nom,2}',['u_1: ' addText],['u_2: ' addText]);
% plot(tVec,distTraj);
% xlabel('Time (s)');
% ylabel('$\|x-x^\star\|/\|x_0-x^\star_0\|$','interpreter','latex');

subplot(2,2,4)
plot(times, energyTraj);
% plot(tVec, energyTraj);
xlabel('Time (s)');
ylabel('Riemann energy')
legend(addText)

figure(2);clf;
hold on;
% if sim_config.include_dist == 1
%     visualize_dist_area(Dist_distribution);
% end
if trajGen_config.include_obs == 1
    visualize_obs(trajGen_config.obs,gray_color);
end
% h2 = plot(xTraj(1,:),xTraj(2,:),'r-','Linewidth',1.5);
h2 = plot(xTraj(1,:),xTraj(2,:),'r-','Linewidth',1.5);
h1 = plot(xnomTraj(1,:),xnomTraj(2,:),'k--','Linewidth',1.5);

scatter(xnomTraj(1,1),xnomTraj(2,1),'ko')
scatter(xnomTraj(1,end),xnomTraj(2,end),'k*')
scatter(xTraj(1,end),xTraj(2,end),'r*')
axis square
xlabel('p_x (m)')
ylabel('p_z (m)')
legend([h1,h2],{'Planned','Actual'});
xlim([0 8]);
ylim([0 8]);
w_max = 1;
if sim_config.save_sim_rst == 1
%     file_name = ['newmodel_adaptiveCCM_ode1_constantintensity_015_03_00_88_gain10'];
    file_name = ['backup2_adaptiveCCM_ode23_constantintensity_03_00_88_gain100'];
    save(file_name,'times','xTraj','uTraj','xnomTraj','unomTraj','energyTraj','dist_config','sim_config','plant','controller');
end