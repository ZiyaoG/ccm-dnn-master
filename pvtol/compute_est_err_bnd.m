% compute the Lipschitz bound of the disturbance
norm_jac = nan*ones(1,1e6);
ii = 1;
for px = 0:0.25:10
    for pz=0:0.25:10
        for vx = -2:0.1:2
            for vz = -2:0.1:2
                x = [px pz 0 vx vz 0]';
                [force,jac] = distLearned(x);
                norm_jac(ii) = norm(jac);
                ii = ii+1;
            end
        end
    end
end
norm_jac = norm_jac(~isnan(norm_jac));
max(norm_jac)

Ld = 5; bd = 0;
theta = Ld*norm([15 15 pi/3 2 1 pi/3])+bd;
max_Bd = norm(plant.B)*theta;

f_norm_vecs = nan*ones(1,1e6);
ii = 1;
for vx = -2:0.2:2
    for vz = -2:0.2:2
        for phi = -pi/3:0.1:pi/3
            for phidot = -pi/3:0.1:pi/3
                x = [0 0 phi vx vz phidot]';
                fx = plant.f_fcn(x); 
                f_norm_vecs(ii) = norm(fx);
                ii=ii+1;
            end
        end
    end
end
max_f = max(f_norm_vecs(~isnan(f_norm_vecs)));
max_Bu = norm(plant.B)*norm([umax; umax]);
phi = max_f+max_Bu+max_Bd;

alpha_T = 2*sqrt(plant.n)*Ld*phi*distEst_config.Ts+sqrt(plant.n)*(1-exp(-distEst_config.a_pred*distEst_config.Ts))*theta;
delta_T1 = norm(pinv(plant.B))*theta
delta_T2 = norm(pinv(plant.B))*alpha_T