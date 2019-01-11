clear;
close all;
tic;
%% parameters
N=3;
% 1->2, 2->3, 3->3
H_o1 = [-eye(2) eye(2) zeros(2)];
H_o2 = [zeros(2) -eye(2) eye(2)];
H_o3 = [zeros(2) zeros(2) eye(2)];
% after modification
R_rel = blkdiag(0.05^2,0.05^2);
R_abs = blkdiag(0.1^2,0.1^2);
% H_m = [sqrt(R_rel^-1)*H_o1;sqrt(R_rel^-1)*H_o2;sqrt(R_abs^-1)*H_o3];

V_i = 0.5; delta = 0.1;
Q_max_i = [eye(2)*(0.1*V_i)^2];
Q_bar = blkdiag(Q_max_i,Q_max_i,Q_max_i);

%% Solve the MARE by brute force iteration
f_q = @(r) r/(1-r);
f_snr = @(L) blkdiag(f_q(L(1))*eye(2),f_q(L(2))*eye(2),f_q(L(3))*eye(2));
f_W = @(L) ones(6) + f_snr(L)^-1*blkdiag(ones(2),ones(2),ones(2));

load P_CEKF.mat;
pi_init = P_CEKF{floor(50/delta)+1};

f_H = @(F) [sqrt((F(1)*R_rel)^-1)*H_o1;sqrt((F(2)*R_rel)^-1)*H_o2;sqrt((F(3)*R_abs)^-1)*H_o3];
Pi_ss = @(F,L) tr_MARE(pi_init,delta^2*Q_bar,f_H(F),f_W(L));
Pi_ss_opt = @(X) tr_MARE(pi_init,delta^2*Q_bar,f_H(X(1:3)),f_W(X(4:6)));
%% OPT - solver-based
tic;
mu_o = 40; mu_c = 200;
% f_cost = @(F,L) mu_o*sum(F.^-1)/delta + mu_c*max(max(diag(F.^-1/delta)*diag(L)));
f_cost = @(F,L) mu_o*sum(F.^-1)/delta + mu_c*trace(diag(F.^-1/delta)*diag(L));
f_cost_opt = @(X) f_cost(X(1:3),X(4:6));

X_init = [5 5 5 0.5 0.5 0.5]; cost_init = f_cost_opt(X_init);
nonlcon_opt = @(X)deal(Pi_ss_opt(X)-0.9*Pi_ss_opt(X_init),[]);
[X_opt,cost_opt] = ga(f_cost_opt,6,[],[],[],[],[1 1 1 0.01 0.01 0.01],[20 20 20 0.99 0.99 0.99],nonlcon_opt,[1 2 3])
% [X_opt,cost_opt] = fmincon(f_cost_opt,X_init,[],[],[],[],[1 1 1 0.6],[1 1 1 0.99],nonlcon_opt,[1 2 3])
ratio_tr = Pi_ss_opt(X_opt)/Pi_ss_opt(X_init)
ratio_cost = cost_opt/f_cost_opt(X_init)
toc