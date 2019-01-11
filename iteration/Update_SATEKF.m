%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%created by Solmaz S. Kia%%%%%%%%%%%%%%%
%%%%%%%%%%Last Revised April 2014%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%UCSD%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%modified by Qi Yan%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Last Revised August 2018%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r_a_bar_send,Gamma_col_send,DeltaP_send,param] = Update_SATEKF(a,b,X_a,X_b,phi_true,X_hat_a,P_a,phi_a,X_hat_b,P_b,phi_b,k,noises,coef_Z)
%% load dt and N when this function is firstly used
persistent dt N n_x fun_H_a fun_H_b C_mat Gamma_s PI_s gamma_p PI;
if isempty(N)
    [~,~,~,N] = RobotInit();
    n_x=2*N;
    Gamma_s = cell(1,N);
    PI_s = cell(1,N);
    PI = cell(1,N*(N-1)/2);
    gamma_p = 1;
    for i=1:N
        Gamma_s{i} = zeros(2);
        PI_s{i} = zeros(2);
    end
    for i=1:N*(N-1)/2
        PI{i} = zeros(2);
    end
end
if isempty(dt)
    [dt,~,~,~,~,~,~,coef_Z] = IterationInit();
    coef_Z = coef_Z(end);
end
if (isempty(fun_H_a) || isempty(fun_H_b))
    C_mat = @(x)[cos(x) -sin(x);sin(x) cos(x)];
    fun_H_a = @(x)-eye(2)*C_mat(x)';
    fun_H_b = @(x) eye(2)*C_mat(x)';
end
%% collect measurement data
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
    phi_est_a =    phi_true(a);
%     z_ab_est  = [x_est_a;y_est_a;phi_est_a];
    z_ab_est  = [x_est_a;y_est_a];
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
%% Calcualte S_ab
% P_check_a = P_a - (1-gamma_p)*phi_a*Gamma_s{a}*Gamma_s{a}'*phi_a';
% P_check_b = P_b - (1-gamma_p)*phi_b*Gamma_s{b}*Gamma_s{b}'*phi_b';
P_check_a = P_a - (1-gamma_p)*phi_a*PI_s{a}*phi_a';
P_check_b = P_b - (1-gamma_p)*phi_b*PI_s{b}*phi_b';
if gamma_p == 0
end
if a ~= b % inter-robot measurement
    H_a = fun_H_a(phi_true(a));
    H_b = fun_H_b(phi_true(b));
    if a > b
        Pi_ba = PI{find_index(b,a)};
        Pi_ab = Pi_ba';
    elseif a < b
        Pi_ab = PI{find_index(a,b)};
        Pi_ba = Pi_ab';
    end
    R_a = R_a(1:2,1:2);
    S_ab = R_a + H_a*P_check_a*H_a' + H_a*phi_a*Pi_ab*phi_b'*H_b' ...
               + H_b*P_check_b*H_b' + H_b*phi_b*Pi_ba*phi_a'*H_a';
elseif a == b
    H_a = eye(2);
    R_a = R_a(1:2,1:2);
    S_ab = R_a + H_a*P_check_a*H_a';
end
S_ab = S_ab + (1-gamma_k)*coef_Z^-1;
%% Scheduled measurements
% gamma_k = 1; s_gamma = 1;
% % R = 0.15;
% % R = 0.5;
% m = 2; eta = norminv(1-(1-(1-R)^(1/m))/2);
% Sigma_k = chol(S_ab,'lower'); r_k = Sigma_k^-1*r_a;
% if max(abs(r_k)) < eta
%     gamma_k = 0;
%     s_gamma = sqrt(2/pi)*eta*exp(-eta^2/2)/(2*normcdf(eta)-1);
% end
param = [gamma_k,h_k];
%% Continuing calculation at the server
Gamma_col = cell(1,N);
if a ~= b % inter-robot measurement
    for i = 1:N
        if i == a
            Gamma_col{i} = (phi_a^-1*P_check_a*H_a' + Pi_ab*phi_b'*H_b')*S_ab^-0.5;
        elseif i == b
            Gamma_col{i} = (phi_b^-1*P_check_b*H_b' + Pi_ba*phi_a'*H_a')*S_ab^-0.5;
        else % not a or b
            if i > a
                Pi_ia = PI{find_index(i,a)}';
            elseif i < a
                Pi_ia = PI{find_index(i,a)};
            end
            if i > b
                Pi_ib = PI{find_index(i,b)}';
            elseif i < b
                Pi_ib = PI{find_index(i,b)};
            end
            Gamma_col{i} = (Pi_ia*phi_a'*H_a' + Pi_ib*phi_b'*H_b')*S_ab^-0.5;
        end
    end
elseif a == b % absolute positioning
    for i = 1:N
        if i == a
            Gamma_col{i} = (phi_a^-1*P_check_a*H_a')*S_ab^-0.5;
        else % not a or b
            if i > a
                Pi_ia = PI{find_index(i,a)}';
            elseif i < a
                Pi_ia = PI{find_index(i,a)};
            end
            Gamma_col{i} = (Pi_ia*phi_a'*H_a')*S_ab^-0.5;
        end
    end
end
r_a_bar = S_ab^-0.5*r_a;
%% update the variables stored at server
% gamma_k = 1;
if gamma_k == 0 % if this is under threshold update
    r_a_bar_send = zeros(size(r_a_bar));
    gamma_p = 0;
%     Gamma_s = Gamma_col;
    Gamma_col_send = cell(1,N); DeltaP_send = cell(1,N);
    for i=1:N
        Gamma_col_send{i} = zeros(2);
        DeltaP_send{i} = zeros(2);
        PI_s{i} = PI_s{i} + Gamma_col{i}*Gamma_col{i}';
    end
%     for i=1:N
%         Gamma_s{i} = zeros(2);
%     end
%     Gamma_col_send = Gamma_col;
elseif gamma_k == 1 % if this is above thresold update
    gamma_p = 1;
    Gamma_col_send = Gamma_col;
    r_a_bar_send = r_a_bar;
    DeltaP_send = PI_s;
    for i=1:N
        PI_s{i} = zeros(2);
    end
end
for i=1:N
%     PI_s{i} = PI_s{i} - Gamma_col{i}*Gamma_col{i}';
    for j=i+1:N
       PI{find_index(i,j)} = PI{find_index(i,j)} - Gamma_col{i}*Gamma_col{j}';
    end
end

if k == 1500
end
end

function index = find_index(a,b)
persistent N data_index;
if isempty(N)
    [~,~,~,N] = RobotInit();
    data_index = zeros(N);k = 1;
    for i = 1:N
        for j = i+1:N
            data_index(i,j) = k;
            k = k + 1;
        end
    end
end
% find the PI index for correlation between a,b
if a == b
    disp('a = b is not allowed!')
    return;
elseif a>b % so that a < b always holds
    tmp = a; a = b; b = tmp;
end
index = data_index(a,b);
end
