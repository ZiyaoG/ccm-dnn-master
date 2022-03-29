function ua= adaptive(t,xa,plant,controller,dist_config,sim_config)
persistent t_pre_a beq_pre_a copt_pre_a Erem_pre_a thetahat

if isempty(t_pre_a) || (t == 0 && t_pre_a ~=0)
    t_pre_a = -3;
    Erem_pre_a = Inf;
end
if isempty(thetahat)
    thetahat = [0;0];
end
geodesic = controller.geodesic; 
x_nom = controller.x_nom_fcn(t);
u_nom = controller.u_nom_fcn(t);
if isempty(geodesic) 
    % when the CCM metric is constant and the geodesic is a straight line 
%     gamma ;
    gamma = [x_nom xa];
    gamma_s = xa-x_nom;
    Erem = gamma_s'*controller.M*gamma_s; 
    gamma_s1_Mx = gamma_s'*controller.M;
    gamma_s0_Mxnom = gamma_s1_Mx;
else 
    n = plant.n; N = geodesic.N; D = geodesic.D; 
    if isempty(beq_pre_a)
        beq_pre_a = zeros(2*plant.n,1);
        copt_pre_a = zeros(plant.n*(geodesic.D+1),1);
    end
    %------------ for testing -----------------
    % Erem = 0;
    % ue = [u;Erem];
    % return;
    %------------------------------------------

    beq = [x_nom;xa];    
    % get the initial value of c corresponding to a straight line
    c0 = zeros(n*(D+1),1);
    %     for i =1:n    
    %         c0((i-1)*(D+1)+1,1) = xStar(i);
    %         c0((i-1)*(D+1)+2,1) = x(i)- xStar(i);    
    %     end
    % vectorized format to improve computational efficiency
    i =1:n;    
    c0((i-1)*(D+1)+1,1) = x_nom;
    c0((i-1)*(D+1)+2,1) = xa - x_nom; 

    % tic;
    if norm(beq-beq_pre_a)<1e-8 && ~isinf(Erem_pre_a)
        copt = copt_pre_a;
        Erem = Erem_pre_a;
    else
        % ----------------- use OPTI -----------------------------------------
        % Opt = opti('fun',geodesic.costf,'grad',geodesic.grad,'eq',geodesic.Aeq,beq,'ndec',geodesic.ndec,'x0',c0,'options',geodesic.opts_opti);
        % [copt,Erem,exitflag,info] = solve(Opt,c0);
        %     nlprob = convIpopt(Opt.prob,geodesic.opts_opti);    
        %     nlprob = convMatlab(Opt.prob,geodesic.opts_opti); 
        % --------------------------------------------------------------------

        % --------------- ipopt ----------------------------------------------
        % geodesic.nlprob.options.rl = beq;
        % geodesic.nlprob.options.ru = beq;
        % geodesic.nlprob.x0 = c0;
        % [copt,Erem,exitflag,info] = opti_ipopt(geodesic.nlprob,c0);
        % --------------------------------------------------------------

        % ---------------- matlab -----------------------
        geodesic.nlprob.beq = beq;
        geodesic.nlprob.x0 = c0;
        [copt,Erem,exitflag,info] = fmincon(geodesic.nlprob);
        if exitflag<0
            disp('geodesic optimization problem failed!');
        end
        % ------------------------------------------------
        beq_pre_a = beq;
        copt_pre_a = copt;
        Erem_pre_a = Erem;
    end
    % toc;
    % ----------------- compute the control law -----------------------
    %     tic;
    %     gamma = zeros(n,N+1);
    %     gamma_s = zeros(n,N+1);  
    %     for i = 1:n   
    %        gamma(i,:) = copt((i-1)*(D+1)+1:i*(D+1),:)'*T;       % gamma(i) is 1*(N+1); the ith elment of gamma on all the (N+1) nodes
    %        gamma_s(i,:) = copt((i-1)*(D+1)+1:i*(D+1),:)'*T_dot;
    %     end  
    %     toc;
    % vectorized format (more computationally efficient)
    copt = transpose(reshape(copt,D+1,n)); % the ith row corresponds to the ith element
    gamma = copt*geodesic.T;
    gamma_s = copt*geodesic.Tdot;
% ----------------------------------------------------------------

% % -------- verify whether the curve found is really a geodesic ----------
% % according to equation (11) in Leung  & Manchester
% error = 0;
% for k=1:N+1
%     error = error + (gamma_s(:,k)'*(controller.W_fcn(gamma(:,k))\gamma_s(:,k))-Erem)^2*geodesic.w_cheby(k);
% end
% error = sqrt(error)/Erem;
% if error>=1e-5
% %     disp('The curve optimized is probably not a geodesic!');
%     fprintf(1,'t= %.2e, Error = %.3e, the curve optimized is probably not a geodesic!\n',t,error);
%     if error> 1e-2
%         pause;
%     end
% end
% % -----------------------------------------------------------------------
end
% tic;

% --------------------------- adaptive ------------------------------------
Gamma = 1e1;
phi = (xa(4)^2+xa(5)^2)/((xa(1)-dist_config.center(1))^2+(xa(2)-dist_config.center(2))^2+1);
gamma_s1 = gamma_s(:,end);
dthetahat = -Gamma*phi*plant.B_fcn(xa)'/controller.W_fcn(xa)*gamma_s1;
thetahat = thetahat + dthetahat*sim_config.step_size;
uu = u_nom + phi*thetahat;

ua = [uu;Erem];
if (t-t_pre_a>= 0.4) && mod(t,1)< 0.1
    fprintf('t = %.1f s\n',t);
    t_pre_a = t;
end
% toc;
end