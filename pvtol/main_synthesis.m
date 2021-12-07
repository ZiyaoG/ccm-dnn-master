%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Synthesis of a CCM for control of a 3D quadrotor

%  Author: Pan Zhao, UIUC, Advanced Controls Research Lab,
%  panzhao2@illinois.edu
% 
%  Last update: Nov 19, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
%% 
yalmip('clear');
controller.type = CtrlDesignOpts.ccm;       % ensure "ccm" is selected
controller.ccm_law_form = ...
    CtrlDesignOpts.ccm_min_norm;            % ensure that "ccm_min_norm" is selected
controller.ccm_mat_ineq_form = ...
    CtrlDesignOpts.ccm_mat_ineq_use_B_perp; % For ccm approach: ccm_mat_ineq_use_rho, or ccm_mat_ineq_use_B_perp . 
controller.opt_pos_dev = 1;                 % whether to focus on mitigate the effect of disturbance on the position states
lambda = 0.8;                               % a line search or bisection search should be performed for lambda
consider_state_set = 1;                     % whether to consider a compact set for the states when formulating the constraints
W_lower_bound = 1e-2;
biSearch_lambda = 0;                        % 1 for conducting a bisection search for lambda in designing the robust CCM controller 
maxNum_bs = 10;                             % maximum number of trials for bisection-search 
lambdaSearchRange = [0.1 10^1];             % search range for lambda
save_rsts = 1;                              % whether to save the results to a file

if controller.type ~= CtrlDesignOpts.ccm || controller.ccm_law_form ~= CtrlDesignOpts.ccm_min_norm
    error('Please ensure that the controller type is CCM and the control law form is min-norm');    
end
%----------------- settings for searching CCM & RCCM ----------------------
n = 6;nw = 1; nu =2;
x = sdpvar(n,1); x_store = x;
W_states_index = [3 4];
if controller.opt_pos_dev == 0
      C= [eye(n); zeros(nu,n)]; D =  [zeros(n,nu); 1*eye(nu)];  % consider both states and inputs
elseif controller.opt_pos_dev == 1
    C = [eye(2) zeros(2,4);zeros(nu,n)]; 
    D = [zeros(2,nu); 1*eye(nu)];                               % consider only position states and inputs
end
% ------------------ Stateconstraints for metric synthesis ----------------
p_lim = pi/3;  % phi
pd_lim = pi/3; % phi_dot 
vx_lim = 2;    % vx
vz_lim = 1;    % vz
state_set.box_lim = [p_lim^2-x(3)^2; vx_lim^2-x(4)^2; pd_lim^2-x(6)^2;  vz_lim^2-x(5)^2]*0.001;
state_set.num_consts_4_W_states = 2;        % # constraints from box_lim that involve states on which the metric W depends
state_set.other_lim_states = [x(6);x(5)]; 
state_set.lagrange_deg_W = 4;               % for the bounds of W
state_set.lagrange_deg_ccm = 4;             % for ccm condition
state_set.p_lim = p_lim;
state_set.pd_lim = pd_lim;
state_set.vx_lim = vx_lim;
state_set.vz_lim = vz_lim;
% -------------------------------------------------------------------------

% ------------------ load system parameters -------------------------------
load_system_parameters; 

% approximating sin_x/cos_x with Chebshev polynomials
sinx = @(x) 0.9101*(x./(pi/3)) - 0.04466*(4*(x./(pi/3)).^3 - 3*(x./(pi/3))); % 0.8799 (close to 0.9101*3/pi, -0.03915
cosx = @(x) 0.7441 -0.2499*(2*(x./(pi/3)).^2 -1);                            % 0.7652, -0.2299

% ------------------------- system dynamics -------------------------------
sin_p = sinx(x(3)); cos_p = cosx(x(3));
f = [x(4)*cos_p - x(5)*sin_p;    %px
    x(4)*sin_p + x(5)*cos_p;     %pz
    x(6);                        %phi
    x(6)*x(5)-plant.g*sin_p;         %vx
    -x(6)*x(4)-plant.g*cos_p;        %vz
    0];                          %phi_dot
% f written as a function handle: can also work when x has multiple columns
f_fcn = @(x) [x(4,:).*cos(x(3,:)) - x(5,:).*sin(x(3,:));    %px
            x(4,:).*sin(x(3,:)) + x(5,:).*cos(x(3,:));     %pz
            x(6,:);                        %phi
            x(6,:).*x(5,:)-plant.g*sin(x(3,:));         %vx
            -x(6,:).*x(4,:)-plant.g*cos(x(3,:));        %vz
            zeros(1,size(x,2))];          
        
B = [zeros(4,2); 1/plant.m 1/plant.m; plant.l/plant.J -plant.l/plant.J]; 
B_perp = [eye(4); zeros(2,4)];
df_dx = jacobian(f,x);
A = df_dx;
nz = size(C,1);

%%%%%%%%%% for testing with(Chebyshev polynomials) approximated fcns%%%%%%%
% using Chebyshev polynomial approximation.
% approximated f_approx, df_dx_fcn
x = x_store;
s = sdisplay(f);
s2 = sdisplay(df_dx);
syms x [n 1]
syms f_approx_fcn [n 1]
syms df_dx_approx_fcn [n n]
for i=1:n    
    f_approx_fcn(i,1) = eval(s{i});    
    for j=1:n
        df_dx_approx_fcn(i,j) = eval(s2{i,j});
    end
end
f_approx_fcn = matlabFunction(f_approx_fcn,'Vars',{x});
df_dx_approx_fcn = matlabFunction(df_dx_approx_fcn,'Vars',{x});
x = x_store;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plant.sinx = sinx;
plant.cosx = cosx;
plant.df_dx = df_dx;
plant.f_fcn = f_fcn;
plant.A = A;
plant.B = B;
plant.B_fcn = @(x) B;
plant.Bpinv_fcn = @(x) pinv(B);
plant.dynamics = @(x,u) f_fcn(x)+ B*u;

plant.B_perp = B_perp;
plant.C = C;
plant.D = D;
plant.n = n; plant.nu=nu; plant.nz = nz; 

W_states = x(W_states_index);
f_phi_fcn = @(x) x(6);
f_vx_fcn = @(x) x(6)*x(5)-plant.g*sin(x(3));
Bw_phi_fcn = @(x) 0;
Bw_vx_fcn = @(x) cos(x(3));

%%%%%%%% for testing using (Chebyshev polynomials) approximated fcns%%%%%%%
f_vx_approx_fcn = @(x) x(6)*x(5)-plant.g*sinx(x(3));
Bw_vx_approx_fcn = @(x) cosx(x(3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v_W = monolist(W_states,4);         % monomials of phi and vx up to degree 4
dv_W_dx = jacobian(v_W,W_states);
n_monos_W = length(v_W);
W_coef = sdpvar(n,n,n_monos_W);
W = zeros(n);
for i=1:n_monos_W
    W = W+ W_coef(:,:,i)*v_W(i);
end

dv_W_dt = dv_W_dx*f(W_states_index);


dW_dt = zeros(n);
for i=1:n_monos_W
    dW_dt = dW_dt+ W_coef(:,:,i)*dv_W_dt(i); 
end
[rho,c_rho,v_rho] = polynomial(W_states,2);
controller.rho = rho;
controller.c_rho = c_rho; 
controller.v_rho = v_rho;
n_monos_rho = length(v_rho);

controller.W_lower_bound = W_lower_bound;
state_set.consider_state_set = consider_state_set;
state_set.W_states = W_states;
state_set.W_states_index = W_states_index;

%% Search for CCM
% YALMIP (SOS) and Mosek are used to solve the SOS problem
% W0 = sdpvar(n,n);
controller.lambda = lambda;%     
%     residual = res
paras_W = W_coef(:);    
[cond_num_W,w_upper,w_lower,W_bar,max_res] = ccm(plant,controller,W,dW_dt,paras_W,lambda,state_set);


% ---------using a constant W matrix for testing or debugging----------
% ---------e.g., whether the optimized geodesic is a straight line.---- 
% W_fcn = @(x) Pinv; 
% dW_fcn = {@(x) zeros(n), @(x) zeros(n), @(x) zeros(n)};
% ---------------------------------------------------------------------

controller.w_upper = w_upper;
controller.w_lower = w_lower;    
controller.W_bar = W_bar;

% ------------------------- extract W_fcn & dW_fcn ------------------------
W_coef = value(W_coef);  
W_coef(abs(W_coef)<=1e-10) = 0;
x = x_store; % must ensure that v_W and s contain "x" instead of "x_store"
W_fcn = zeros(n);
for i=1:n_monos_W
    W_fcn = W_fcn+ W_coef(:,:,i)*v_W(i);
end
W_fcn = clean(W_fcn, 1e-10);
dv_W_dx = clean(dv_W_dx, 1e-10);
s = sdisplay(W_fcn);
s2 = sdisplay(dv_W_dx);
syms x [n 1]
syms W_fcn [n n]
for i=1:n
    for j=1:n
        W_fcn(i,j) = eval(s{i,j});
    end
end
%     W_fcn = arrayfun(@(si) eval(si), s);
matlabFunction(W_fcn,'File','W_fcn1','Vars',{x});
W_fcn = matlabFunction(W_fcn,'Vars',{x});
% W_fcn = @(x) W_fcn1(x(3),x(4));


[n1,n2]= size(dv_W_dx);
syms dv_W_dx_sym [n1 n2]

for i=1:n1
    for j=1:n2
        dv_W_dx_sym(i,j) = eval(s2{i,j});
    end
end    
dW_dphi = zeros(n);
dW_dvx = zeros(n);
for i=1:n_monos_W
    dW_dphi = dW_dphi + W_coef(:,:,i)*dv_W_dx_sym(i,1); 
    dW_dvx = dW_dvx + W_coef(:,:,i)*dv_W_dx_sym(i,2);
end
%     dW_dx_fcn = (i==3)*dW_fcn1 + (i==4)*dW_fcn2;
%     s = sdisplay(dW_dx_fcn);
matlabFunction(dW_dphi,'File','dW_dphi','Vars',{x});
matlabFunction(dW_dvx,'File','dW_dvx','Vars',{x});
dW_dphi = matlabFunction(dW_dphi,'Vars',{x});
dW_dvx = matlabFunction(dW_dvx,'Vars',{x});
dW_dxi_fcn = @(i,x) (i==3)*dW_dphi(x)+(i==4)*dW_dvx(x);


dW_dt_fcn = @(x) dW_dphi(x)*f_phi_fcn(x) + dW_dvx(x)*f_vx_fcn(x); 

controller.W_fcn = W_fcn;
controller.dW_dxi_fcn = dW_dxi_fcn;
controller.dW_dt_fcn = dW_dt_fcn;

%% --------------check CCM conditions and save data----------
lambda0 = lambda;
if isfield(plant,'df_dx')
    plant= rmfield(plant,{'df_dx','A'});
end
if isfield(controller,'rho')
    controller= rmfield(controller,{'rho','c_rho','v_rho'});
end
if isfield(state_set,'box_lim')
    state_set = rmfield(state_set,{'box_lim','other_lim_states','W_states'});
end
plant.state_set = state_set;

disp('Checking CCM conditions ...');

lambda = 0.998*lambda; % slightly reduce lambda to avoid infeasible problems due to numerical issues    
ctrl_N = 10;
p_range = linspace(-p_lim, p_lim, ctrl_N);
vx_range = linspace(-vx_lim, vx_lim, ctrl_N);
vz_range = linspace(-vz_lim, vz_lim, 8);
pd_range = linspace(-pd_lim, pd_lim, 8);

df_dx_fcn = @(x) [0,0,-x(4)*sin(x(3))-x(5)*cos(x(3)),cos(x(3)),-sin(x(3)),0; 
               0,0, x(4)*cos(x(3))-x(5)*sin(x(3)),sin(x(3)), cos(x(3)),0;
               zeros(1,5),1;
               0,0,-plant.g*cos(x(3)),0,x(6),x(5);
               0,0, plant.g*sin(x(3)),-x(6),0,-x(4);
               zeros(1,6)];
A_fcn = @(x,w)  df_dx_fcn(x) + dBw_dx_fcn(x)*w;   
A_approx_fcn = @(x,w)  df_dx_approx_fcn(x) + dBw_dx_approx_fcn(x)*w;   
controller.df_dx_fcn = df_dx_fcn;
controller.A_fcn = A_fcn;
controller.A_approx_fcn = A_approx_fcn;
 
% --------------- following the approach in the following paper -------
% S.  Singh, et al. Robust  feedback  motion  planning  via
% contraction  theory. IJRR, 2019.
delta_u = zeros(ctrl_N,ctrl_N,ctrl_N,ctrl_N);
eig_CCM = zeros(ctrl_N, ctrl_N, ctrl_N, ctrl_N);
eig_W = zeros(ctrl_N,ctrl_N,2);

for i = 1:length(p_range)
    for j = 1:length(vx_range)
        x = [randn(2,1);p_range(i);vx_range(j);0;0];                
        W = W_fcn(x); % note that W depends only on phi and vx
        M = W\eye(n);
        eig_W(i,j,1) = min(eig(W));
        eig_W(i,j,2) = max(eig(W));
        for k = 1:length(vz_range)
            for l = 1:length(pd_range)
                x = [randn(2,1);p_range(i);vx_range(j);vz_range(k);pd_range(l)];                
                Theta = chol(M);
                L = chol(W);                
%                     f = f_part_fcn(x);

                df_dx0 = df_dx_fcn(x);                
                F = dW_dt_fcn(x) - df_dx0*W - (df_dx0*W)'-2*lambda*W;

                L_inv = inv(L);                
                delta_u_den = eig((L_inv)'*(B*B')*L_inv);
                delta_u(i,j,k,l) = 0.5*max(eig((L_inv)'*F*L_inv))/...
                    sqrt(min(delta_u_den(delta_u_den>0)));

                R_CCM = B_perp'*F*B_perp;
                eig_CCM(i,j,k,l) = min(eig(R_CCM));
            end
        end
    end
 
end
Wbar_inv = eye(n)/W_bar;
Wbar_inv_eigs = eig(Wbar_inv(1:2,1:2));

fprintf(1,'CCM, lambda = %.2f, cond(W) = %.3f\n',lambda,cond_num_W);
fprintf(1,'min and max eigenvalues of W: %.3e, %.3e\n', min(min(eig_W(:,:,1))),max(max(eig_W(:,:,2))));
fprintf(1,'minimum eigenvalue of CCM matrix (should be positive): %.3e\n',min(eig_CCM(:)));

if controller.ccm_mat_ineq_form == CtrlDesignOpts.ccm_mat_ineq_use_rho       
    % extract rho_fcn
    c_rho = value(c_rho);     
    x = x_store; % must ensure that v_W and s contain "x" instead of "x_store"
    rho_fcn = 0;
    for i=1:n_monos_rho
        rho_fcn = rho_fcn+ c_rho(i)*v_rho(i);
    end
    rho_fcn = clean(rho_fcn, 1e-8);
    s = sdisplay(rho_fcn);
    syms x [n 1]
    syms rho_fcn
    rho_fcn = eval(s{1});
    rho_fcn = matlabFunction(rho_fcn,'Vars',{x});
    controller.rho_fcn = rho_fcn;
end
if save_rsts == 1
    file_name = ['ccm_' num2str(lambda0) '_plim_' num2str(p_lim/pi,2) 'pi.mat'];
    save(file_name,'plant','controller','state_set');
end
