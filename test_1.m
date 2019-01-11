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
H_bar = [H_o1;H_o2;H_o3];
R_rel = blkdiag(0.05^2,0.05^2);
R_abs = blkdiag(0.1^2,0.1^2);
R_bar = blkdiag(R_rel,R_rel,R_abs);

V_i = 0.5; delta = 0.1;
Q_max_i = [eye(2)*(0.1*V_i)^2];
Q = blkdiag(Q_max_i,Q_max_i,Q_max_i);

%%
Pi_ss = @(W)dare(eye(2*N),H_bar',delta^2*Q,W,zeros(2*N,2*N),eye(2*N));
mat_f = @(f) blkdiag(eye(2)/f(1),eye(2)/f(2),eye(2)/f(3));
Pi_bar = @(f,Z)Pi_ss(R_bar*mat_f(f) + Z^-1);

f_rate = @(f,Z) 1 - 1/sqrt(det(eye(2*N)+(H_bar*Pi_bar(f,Z)*H_bar' + R_bar*mat_f(f))*Z));
R_star = @(f,Z) (f_rate(f,Z)*(R_bar*mat_f(f))^-1 + (1-f_rate(f,Z))*(R_bar*mat_f(f) + Z^-1)^-1)^-1;
Pi_exp = @(f,Z) Pi_ss(R_star(f,Z));

% coef_Z = [0.1:0.1:10,10:100];
coef_Z = [1:100];
coef_f_o = 0.1:0.1:1;
TR_bar = zeros(10,10); TR_exp = zeros(10,10); Rate = zeros(10,10);
Cost = zeros(10,10);
mu_o = 50; mu_c = 50;
for j = 1:length(coef_f_o)
    f = ones(1,3)*coef_f_o(j);
    for i = 1:length(coef_Z)
        Z = eye(2*N)*coef_Z(i);
        TR_bar(i,j) = trace(Pi_bar(f,Z));
        TR_exp(i,j) = trace(Pi_exp(f,Z));
        Rate(i,j) = f_rate(f,Z);
        Cost(i,j) = mu_o*f(1) + mu_c*f(1)*Rate(i,j);
    end
end
%% OPT - solver-based
% f_cost = @(X) mu_o*X(1) + mu_c*X(1)*f_rate(X(1),X(2)*eye(2*N)); % X(1) - rate; X(2) - Z
mu_o = 40; mu_c = 150;
f_cost = @(X) mu_o*sum(X(1:3)) + mu_c*f_rate(X(1:3),X(4)*eye(2*N))*max(X(1:3)); % X(1:3) - observation rate; X(4) - z
X_init = [0.5 0.5 0.5 40]; cost_init = f_cost(X_init);
% [X_opt,cost_opt] = fmincon(f_cost,X_init,[],[],[],[],[0.01 0.01 0.01 0],[1 1 1 200],@nonlcon)
f_cost_int = @(X) f_cost([1/X(1) 1/X(2) 1/X(3) X(4)]);
[X_opt,cost_opt] = ga(f_cost_int,4,[],[],[],[],[1 1 1 1],[10 10 10 100],@nonlcon_int,[1 2 3])
f_fate_int = @(X) f_rate([1/X(1) 1/X(2) 1/X(3)],X(4)*eye(2*N));
rate_cur = f_fate_int(X_opt)
tr_cur = trace(Pi_bar(X_opt(1:3),X_opt(4)*eye(2*N)));
%% OPT - problem-based

%% iteration to verify the OPT
pi_cell = cell(10,1); Tr_pi = zeros(10,1);
load P_CEKF.mat;
pi_cell{1} = P_CEKF{floor(50/delta)+1};
f_o = min(X_opt(1:3)); % frequency of observation in a propagation time step
step_o = f_o;
for i=1:1e3
    if mod(i,step_o) == 0
        is_update = 1;
    else
        is_update = 0;
    end
    pi_cell{i+1} = pi_cell{i} - is_update*pi_cell{i}*H_bar'*(H_bar*pi_cell{i}*H_bar' + R_bar)^-1*H_bar*pi_cell{i};
    pi_cell{i+1} = pi_cell{i+1} + delta^2*Q;
    Tr_pi(i) = trace(pi_cell{i});
end
[Pi_best,~,~] = dare(eye(2*N),H_bar',delta^2*Q,R_bar,zeros(2*N,2*N),eye(2*N));
[Pi_low,~,~] = dare(eye(2*N),H_bar',delta^2*Q,R_bar/f_o,zeros(2*N,2*N),eye(2*N));

%%
figure;
hold on;
box on;
plot(Tr_pi);        
line(xlim,[trace(Pi_best) trace(Pi_best)],'linestyle','--','linewidth',2);
line(xlim,[trace(Pi_low) trace(Pi_low)],'linestyle','-','linewidth',2);
legend('Iterative','Analytical-best','Analytical-low');
xlim([0 50]);

%%
figure;
box on;
[mesh_x,mesh_z] = meshgrid(coef_Z,coef_f_o);
% surf(mesh_x,mesh_z,TR_bar');
surf(mesh_x,mesh_z,TR_exp');
ylabel('Frequency');
xlabel('Trace of Z');
zlabel('Uncertainty');

figure;
surf(mesh_x,mesh_z,Rate');
ylabel('Frequency');
xlabel('Trace of Z');
zlabel('Communication rate');

%%
figure; hold on;
box on;
yyaxis left;
% plot(Rate(:,3),TR_bar(:,3),'b-','linewidth',2);
% plot(Rate(:,3),TR_exp(:,3),'b--','linewidth',2);
% plot(Rate(:,8),TR_bar(:,8),'g-','linewidth',2);
% plot(Rate(:,8),TR_exp(:,8),'g--','linewidth',2);
plot(coef_Z,TR_bar(:,3),'b-','linewidth',2);
plot(coef_Z,TR_exp(:,3),'b--','linewidth',2);
plot(coef_Z,TR_bar(:,10),'g-','linewidth',2);
plot(coef_Z,TR_exp(:,10),'g--','linewidth',2);
line(xlim,[trace(Pi_best) trace(Pi_best)],'linestyle','--','linewidth',2);
xlabel('Coefficent of Z');
ylabel('Trace');

yyaxis right;
plot(coef_Z,Rate(:,3),'linewidth',2);
plot(coef_Z,Rate(:,10),'linewidth',2);
ylabel('Com. rate');
legend('Upper bound,f = 0.3','Empirical, f = 0.3','Upper bound,f = 1','Empirical, f = 1','Best bound','Com. rate, f = 0.3','Com. rate, f = 1','location','best');

%%
% figure;
% hold on;
% yyaxis left;
% plot(coef_Z,TR_bar,'linewidth',2);
% plot(coef_Z,TR_exp,'linewidth',2);
% xlabel('Coefficient (trace of Z)');
% ylabel('Trace');
% line(xlim,[trace(Pi_low) trace(Pi_low)],'color','black','linestyle','-.','linewidth',2);
% ylim([0 0.1]);
% xlim([0 100]);
% 
% plot(coef_Z,Rate,'linewidth',2);
% % xlabel('Coefficient (trace of Y)');
% ylabel('Communication rate');
% legend('Upper bound'
% yyaxis right,'A-Lower bound','Lower bound','Rate-lower bound','Rate-upper bound');
toc;
%%
function [c,ceq] = nonlcon(X)
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
% c = trace(Pi_exp(X(1),X(2)*eye(2*N))) - 1.01*trace(Pi_best);
% c = trace(Pi_bar(X(1),X(2)*eye(2*N))) - 1.2*trace(Pi_best);
X_init = [0.3 0.3 0.3 40];
% c = trace(Pi_bar(X(1:3),X(4)*eye(2*N))) - 1*trace(Pi_bar(X_init(1:3),X_init(4)*eye(2*N)));
c = trace(Pi_exp(X(1:3),X(4)*eye(2*N))) - 1*trace(Pi_exp(X_init(1:3),X_init(4)*eye(2*N)));
end