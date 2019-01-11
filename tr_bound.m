function tr_bound = tr_bound(P_init,V_i,delta)
N=3;
C_mat = @(x)[cos(x) -sin(x);sin(x) cos(x)];
H_1 = C_mat(0)'*[-eye(2) eye(2) zeros(2)];
H_2 = C_mat(0)'*[zeros(2) -eye(2) eye(2)];
H_3 = C_mat(0)'*[zeros(2) zeros(2) eye(2)];
% after modification
H_bar = [H_1;H_2;H_3];
R_rel = blkdiag(0.1^2,0.1^2);
R_abs = blkdiag(0.1^2,0.1^2);
R_bar = blkdiag(R_rel,R_rel,R_abs);

% % V_i = 0.3; delta = 0.5;
Q_max_i = [eye(2)*(0.35*V_i)^2];
Q = blkdiag(Q_max_i,Q_max_i,Q_max_i);
pi_cell = cell(10,1); Tr_pi = zeros(10,1);
pi_cell{1} = P_init;
% ci_alpha = 0.98;
% F_hat = 1/sqrt(ci_alpha)*eye(2*N);
% R_hat = blkdiag(0.5*R_bar,ci_alpha/(1-ci_alpha)*R_bar);
for i=1:1e3
    %     pi_cell{i+1} = F_hat*pi_cell{i}*F_hat' + delta^2*Q_max ...
    %                    - F_hat*pi_cell{i}*H_hat'*(R_hat + H_hat*pi_cell{i}*H_hat')^-1*H_hat*pi_cell{i}*F_hat';
        pi_cell{i+1} = (pi_cell{i}^-1 + H_bar'*R_bar^-1*H_bar)^-1;
%     pi_cell{i+1} = pi_cell{i} - pi_cell{i}*H_bar'*(H_bar*pi_cell{i}*H_bar' + R_bar)^-1*H_bar*pi_cell{i};
%     pi_cell{i+1} = pi_cell{i}*(eye(2*N) + H_bar'*R_bar^-1*H_bar*pi_cell{i})^-1;
    pi_cell{i+1} = pi_cell{i+1} + delta^2*Q;
    Tr_pi(i) = trace(pi_cell{i});
end
tr_bound = Tr_pi(end);
% [X_dare,L_dare,G_dare] = dare(F_hat',H_hat',delta^2*Q_max,R_hat,zeros(2*N,2*N),eye(2*N));
% figure;
% plot(Tr_pi);
end