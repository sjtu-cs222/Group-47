%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%created by Solmaz S. Kia%%%%%%%%%%%%%%%
%%%%%%%%%%Last Revised April 2014%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%UCSD%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%modified by Qi Yan%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Last Revised August 2018%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X_Prop,P_Prop,Phi_Prop_SATEKF,G_Prop] = Prop_SAEKF(X_hat_k,phi_mea,P_k_single,Phi_k,k,sigma_V,sigma_phi)
%% load dt and N when this function is firstly used
persistent dt N fun_F fun_G ;
if isempty(N)
%     [~,~,~,N] = RobotInit();
    N = length(X_hat_k);
end
if isempty (dt)
    [dt,~] = IterationInit();
end
if isempty(fun_F)
%     fun_F = @(x,u) [1 0 -u(1)*sin(x(3))*dt; 0 1 u(1)*cos(x(3))*dt; 0 0 1]; % Jacobian matrix for system function (states)
%     fun_G = @(x,u) [dt*cos(x(3)) 0;dt*sin(x(3)) 0;0 dt]; % Jacobian matrix for system function (inputs)
    fun_F = eye(2);
    fun_G = @(v,phi)[cos(phi) -v*sin(phi);sin(phi) v*cos(phi)];
end
%%
[V,W] = TrajGen(k);
Q = cell(1,N);
for i = 1:N
%     Q{i} = diag([sigma_V(i)^2,sigma_W(i)^2]);
    Q{i} = diag([sigma_V(i)^2,(sigma_phi(i)*V(i))^2]);
%     Q{i} = diag(sigma_V(i)^2);
end

%% EKF Propagation robot-wisely
X_Prop = cell(1,N); this_F = cell(1,N); this_G = cell(1,N);
P_single = cell(1,N); Phi_Prop_SATEKF = cell(1,N);
G_single = cell(1,N);
for i = 1:N
    % after going through the motion equation
    X_Prop{i} = X_hat_k{i} + dt*[V(i)*cos(phi_mea(i));V(i)*sin(phi_mea(i))];
    this_F{i} = fun_F;
    this_G{i} = fun_G(V(i),phi_mea(i));
    P_single{i} = this_F{i}*P_k_single{i}*this_F{i}' +  dt^2*this_G{i}*Q{i}*this_G{i}';
    G_single{i} = dt^2*this_G{i}*Q{i}*this_G{i}';
    Phi_Prop_SATEKF{i} = this_F{i}*Phi_k{i};
end
P_Prop = P_single;
G_Prop = blkdiag(G_single{:});
end