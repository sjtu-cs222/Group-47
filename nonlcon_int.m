function [c,ceq] = nonlcon_int(X)
%% parameters
N=3;
% 1->2, 2->3, 3->3
H_o1 = [-eye(2) eye(2) zeros(2)];
H_o2 = [zeros(2) -eye(2) eye(2)];
H_o3 = [zeros(2) zeros(2) eye(2)];
% after modification
H_bar = [H_o1;H_o2;H_o3];
R_rel = blkdiag(0.05^2,0.05^2);
R_abs = blkdiag(0.1^2,0.1^2);
R_bar = blkdiag(R_rel,R_rel,R_abs);

V_i = 0.3; delta = 0.2;
Q_max_i = [eye(2)*(0.1*V_i)^2];
Q = blkdiag(Q_max_i,Q_max_i,Q_max_i);

%%
Pi_ss = @(W)dare(eye(2*N),H_bar',delta^2*Q,W,zeros(2*N,2*N),eye(2*N));
mat_f = @(f) blkdiag(eye(2)/f(1),eye(2)/f(2),eye(2)/f(3));
Pi_bar = @(f,Z)Pi_ss(R_bar*mat_f(f) + Z^-1);

f_rate = @(f,Z) 1 - 1/sqrt(det(eye(2*N)+(H_bar*Pi_bar(f,Z)*H_bar' + R_bar*mat_f(f))*Z));
R_star = @(f,Z) (f_rate(f,Z)*(R_bar*mat_f(f))^-1 + (1-f_rate(f,Z))*(R_bar*mat_f(f) + Z^-1)^-1)^-1;
Pi_exp = @(f,Z) Pi_ss(R_star(f,Z));

%%
ceq = [];
X_init = [0.5 0.5 0.5 40];
X = [1/X(1) 1/X(2) 1/X(3) X(4)];
c = trace(Pi_exp(X(1:3),X(4)*eye(2*N))) - 1*trace(Pi_exp(X_init(1:3),X_init(4)*eye(2*N)));
end