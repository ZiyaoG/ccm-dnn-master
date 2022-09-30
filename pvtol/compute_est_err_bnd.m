% compute the Lipschitz bound of the disturbance
% distLearned=@(x) actual_dist_fcn(x,dist_config.center,dist_config.radius);
% distLearned = @(x) learned_dist_fcn(x,actual_dist_fcn);
norm_jac = nan*ones(1,1e6);
ii = 1;
for px = 0:0.25:10
    for pz=0:0.25:10
        for vx = -2:0.1:2
            for vz = -2:0.1:2
                x = [px pz 0 vx vz 0]';
                [force,jac] = learned_dist_fcn(x,@actual_dist_fcn);
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

function dist_force = actual_dist_fcn(x)
max_damping = 0.5;
center = [5,5]';
radius = 5;

[dist_intensity,~] = dist_distribution(x(1),x(2),center,radius);
% dist_intensity = 0.3;
dist_force_max = (x(3)^2+x(4)^2)*max_damping;
dist_force = [-1; -1]*(dist_intensity.*dist_force_max); 
end

function [intensity,distance_to_center] = dist_distribution(X,Z,center,radius)
distance_to_center = sqrt((X-center(1)).^2 + (Z-center(2)).^2);
intensity = 1./(distance_to_center.^2+1);

% --------------- using an inverse function 2 ----------------
% intensity = 1./(sqrt(distance_to_center.^2)+1);
end


function [dist_force, d_force_dx]= learned_dist_fcn(x,actual_dist_fcn)
[n,N]  = size(x);

d_force_dx = zeros(2,n,N);
% dist_force = zeros(2,N);

x_reduced = x([1 2 4 5],:)';
dist_force = actual_dist_fcn(x_reduced);
if nargout>1
    [n,N] = size(x);
%     tic;
%     for i = 1:N        
%         % ---------------- using an inverse function -----------------
% %         xx = x(:,i);
%         xx_reduced = x_reduced(i,:);
% 
% %       tic;
%         jac = jacobian(prd_dist,xx_reduced);
% %         toc;
%         d_force_dx(:,[1 2 4 5],i) = jac;
%     end
%     toc;

    jac = jacobianN(actual_dist_fcn,x_reduced);
    d_force_dx(:,[1 2 4 5],:) = jac;
end

end