%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%created by Solmaz S. Kia%%%%%%%%%%%%%%%
%%%%%%%%%%Last Revised April 2014%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%UCSD%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%modified by Qi Yan%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Last Revised August 2018%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X_Prop,P_Prop] = Prop_DEKF(X_hat_k,P_single_k,k,sigma_V,sigma_W)
%% load dt and N when this function is firstly used
persistent dt N fun_F fun_G;
if isempty(N)
%     [~,~,~,N] = RobotInit();
    N = length(X_hat_k);
end
if isempty (dt)
    [dt,~] = IterationInit();
end
if isempty(fun_F)
    fun_F = @(x,u) [1 0 -u(1)*sin(x(3))*dt; 0 1 u(1)*cos(x(3))*dt; 0 0 1]; % Jacobian matrix for system function (states)
    fun_G = @(x,u) [dt*cos(x(3)) 0;dt*sin(x(3)) 0;0 dt]; % Jacobian matrix for system function (inputs)
end
%%
[V,W] = TrajGen(k);
Q = cell(1,N);
for i = 1:N
    Q{i} = diag([sigma_V(i)^2,sigma_W(i)^2]);
end

%% EKF Propagation robot-wisely
X_Prop = cell(1,N); P_Prop = cell(1,N);
for i = 1:N
    % after going through the motion equation
    phi = X_hat_k{i}(3);
    X_Prop{i} = X_hat_k{i}...
              + [V(i)*cos(phi)*dt;
                 V(i)*sin(phi)*dt;
                 W(i)*dt];
    
    this_F = fun_F(X_hat_k{i},[V(i),W(i)]);
    this_G = fun_G(X_hat_k{i},[V(i),W(i)]);
    P_Prop{i} = this_F*P_single_k{i}*this_F' + this_G*Q{i}*this_G';
end

end