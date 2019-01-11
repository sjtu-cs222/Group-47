%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%created by Solmaz S. Kia%%%%%%%%%%%%%%%
%%%%%%%%%%Last Revised April 2014%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%UCSD%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%modified by Qi Yan%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Last Revised August 2018%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X_trueProp,phi_true,phi_mea,sigma_V,sigma_phi] = TrueRobotModel(X_k,phi_true,phi_mea,k,noise_v,noise_w)
%% load dt and N when this function is firstly used
persistent dt N
if isempty (dt)
    [dt,~] = IterationInit();
end
if isempty(N)
%     [~,~,~,N] = RobotInit();
    N = length(X_k);
end
%% computation for true robot motion
[V,W] = TrajGen(k);
[sigma_V,sigma_phi] = ProprioVar(V,k);
X_trueProp = cell(1,N);
for i=1:N
    V_true = V(i) - sigma_V(i)*noise_v(i);
    X_trueProp{1,i} = X_k{1,i} + dt*[V_true*cos(phi_true(i)); V_true*sin(phi_true(i))];
    phi_mea(i) = phi_mea(i) + W(i)*dt;
    phi_true(i) = phi_mea(i) + W(i)*dt - sigma_phi(i)*noise_w(i);
end

end
