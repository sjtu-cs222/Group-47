function rate = fun_rate(f,z)
%% parameters
N=3;
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
f_rate = @(f,Z) 1 - 1/sqrt(det(eye(2*N)+(H_bar*Pi_bar(Z,f)*H_bar' + R_bar/f)*Z));
rate = f_rate(f,z*eye(2*N));
end