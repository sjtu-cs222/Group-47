%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%created by Solmaz S. Kia%%%%%%%%%%%%%%%
%%%%%%%%%%Last Revised April 2014%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%UCSD%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%modified by Qi Yan%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Last Revised August 2018%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [K,S_ab,r_a,param] = Update_TEKF(a,b,X_a,X_b,phi_true,X_hat_a,X_hat_b,P_k,k,noises,coef_Z)
%% load dt and N when this function is firstly used
persistent dt N n_x fun_H_a fun_H_b C_mat gamma_p;
if isempty(N)
    [~,~,~,N] = RobotInit();
    n_x=2*N;
end
if isempty(dt)
    [dt,~,~,~,~,~,~,coef_Z] = IterationInit();
    coef_Z = coef_Z(end);
    gamma_p = 1;
end
if (isempty(fun_H_a) || isempty(fun_H_b))
    C_mat = @(x)[cos(x) -sin(x);sin(x) cos(x)];
    fun_H_a = @(x)-eye(2)*C_mat(x)';
    fun_H_b = @(x) eye(2)*C_mat(x)';
end
%%
[R_rel,R_abs] = ExteroVar(k); % measurement noise variance

x_true_a=X_a(1);
y_true_a=X_a(2);
% phi_true_a=0;
phi_true_a=phi_true(a);

x_true_b=X_b(1);
y_true_b=X_b(2);
% phi_true_b=0;
phi_true_b=phi_true(b);

if a~=b
    R_a=R_rel{1,a};
%     v_z_a=[sqrt(R_a(1,1))*noises(1);sqrt(R_a(2,2))*noises(2);sqrt(R_a(3,3))*noises(3)];
    v_z_a=[sqrt(R_a(1,1))*noises(1);sqrt(R_a(2,2))*noises(2)];
elseif a == b
    clear R_a v_z_a;
    R_a=R_abs{1,a}; 
    R_a = blkdiag(R_a,(2/180*57.3)^2);
    v_z_a=[sqrt(R_a(1,1))*noises(4);sqrt(R_a(2,2))*noises(5)];
%     v_z_a=[v_z_a;sqrt(R_a(3,3))*noises(6)];
    v_z_a=[v_z_a];
end

if a~=b
    C = [cos(phi_true_a)  -sin(phi_true_a)
        sin(phi_true_a)   cos(phi_true_a)];    
%     z_ab=[C'*([x_true_b-x_true_a;y_true_b-y_true_a]); phi_true_b-phi_true_a] + v_z_a;
    z_ab=[C'*([x_true_b-x_true_a;y_true_b-y_true_a])] + v_z_a;

    x_est_a   =    X_hat_a(1);%   X_hat_a(1);
    y_est_a   =    X_hat_a(2);%   X_hat_a(2);
%     phi_est_a=0;
    phi_est_a =    phi_true(a);%   X_hat_a(3);
    
    x_est_b   =    X_hat_b(1);%   X_hat_b(1);
    y_est_b   =    X_hat_b(2);%   X_hat_b(2);
%     phi_est_b = 0;
    phi_est_b =    phi_true(b);%   X_hat_b(3);
    C_est = [cos(phi_est_a)  -sin(phi_est_a)
        sin(phi_est_a)   cos(phi_est_a)];
    
%     z_ab_est = [C_est'*([x_est_b-x_est_a;y_est_b-y_est_a]); phi_est_b-phi_est_a];
    z_ab_est = [C_est'*([x_est_b-x_est_a;y_est_b-y_est_a])];
    
elseif a == b
%     z_ab = [x_true_a;y_true_a;phi_true_a] + v_z_a;
    z_ab = [x_true_a;y_true_a] + v_z_a;
    
    x_est_a   =    X_hat_a(1);
    y_est_a   =    X_hat_a(2);
%     phi_est_a = 0;
    phi_est_a =    phi_true(a);
%     z_ab_est = [x_est_a;y_est_a;phi_est_a];
    z_ab_est = [x_est_a;y_est_a];
end

r_a = z_ab - z_ab_est;
r_a = r_a(1:2);

coef_Z = coef_Z*eye(2);
zeta_k = rand;
phi_k = exp(-0.5*r_a'*coef_Z*r_a);
if zeta_k <= phi_k % not send
    gamma_k = 0;
else
    gamma_k = 1;
end
h_k = 1;

if a~=b % inter-robot measurement
    H_a = fun_H_a(phi_true(a));
    H_b = fun_H_b(phi_true(b));
    P_a = P_k(2*a-1:2*a,2*a-1:2*a); P_b = P_k(2*b-1:2*b,2*b-1:2*b);
    P_ab = P_k(2*a-1:2*a,2*b-1:2*b); P_ba = P_k(2*b-1:2*b,2*a-1:2*a);
    R_a = R_a(1:2,1:2);
    S_ab = R_a + H_a*P_a*H_a' + H_b*P_b*H_b' + H_a*P_ab*H_b' + H_b*P_ba*H_a';
    S_ab = S_ab + (1-gamma_k)*coef_Z^-1;
    K = cell(1,N);
    for j = 1:N
        P_ja = P_k(2*j-1:2*j,2*a-1:2*a);
        P_jb = P_k(2*j-1:2*j,2*b-1:2*b); 
        K{j} = (P_ja*H_a' + P_jb*H_b')*S_ab^-1;
    end
    K = cell2mat(K');
elseif a == b
%     H_a = fun_H_a(X_hat_a,X_hat_a);
    H_a = eye(2);
    P_a = P_k(2*a-1:2*a,2*a-1:2*a);
    R_a = R_a(1:2,1:2);
    S_ab = R_a + H_a*P_a*H_a';
    S_ab = S_ab + (1-gamma_k)*coef_Z^-1;
    K = cell(1,N);
    for j = 1:N
        P_ja = P_k(2*j-1:2*j,2*a-1:2*a);
        K{j} = (P_ja*H_a')*S_ab^-1;
    end
    K = cell2mat(K');
end
%% Scheduled measurements
if gamma_p == 0
end
% gamma_k = 1; h_k = 1;
% m = 2; eta = norminv(1-(1-(1-R)^(1/m))/2);
% Sigma_k = chol(S_ab,'lower');
% r_k = Sigma_k^-1*r_a;
% if max(abs(r_k)) < eta
%     gamma_k = 0;
%     h_k = sqrt(2/pi)*eta*exp(-eta^2/2)/(2*normcdf(eta)-1);
% end

param = [gamma_k,h_k];
gamma_p = gamma_k;
end

