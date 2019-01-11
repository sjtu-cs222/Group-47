clear
close all

% Construct observation matrix
f_mat_C = @(phi) [cos(phi) -sin(phi);sin(phi) cos(phi)];
H = cell(7,1);
% H{1} = f_mat_C(pi/3)*[-eye(2) eye(2) zeros(2)];       % 1->2
% H{2} = f_mat_C(pi/3)*[-eye(2) zeros(2) eye(2)];       % 1->3
% H{3} = f_mat_C(-pi/3)*[eye(2) -eye(2) zeros(2)];       % 2->1
% H{4} = f_mat_C(-pi/3)*[zeros(2) -eye(2) eye(2)];       % 2->3
% H{5} = f_mat_C(pi/6)*[zeros(2) eye(2) -eye(2)];       % 3->2
% H{6} = f_mat_C(pi/6)*[eye(2) zeros(2) -eye(2)];       % 3->1
% H{7} = f_mat_C(pi/6)*[zeros(2) zeros(2) -eye(2)];     % 3->3

H{1} = f_mat_C(0)*[-eye(2) eye(2) zeros(2)];       % 1->2
H{2} = f_mat_C(0)*[-eye(2) zeros(2) eye(2)];       % 1->3
H{3} = f_mat_C(0)*[eye(2) -eye(2) zeros(2)];       % 2->1
H{4} = f_mat_C(0)*[zeros(2) -eye(2) eye(2)];       % 2->3
H{5} = f_mat_C(0)*[zeros(2) eye(2) -eye(2)];       % 3->2
H{6} = f_mat_C(0)*[eye(2) zeros(2) -eye(2)];       % 3->1
H{7} = f_mat_C(0)*[zeros(2) zeros(2) -eye(2)];     % 3->3
H_col = vertcat(H{:});
len_m = length(H_col);
num_sensor = len_m/2;

% System matrix
mat_A = eye(6);
% system process noise
V_i = 20; delta = 0.1;
Q_max_i = eye(2)*(0.35*V_i)^2*delta^2;
Q_bar = blkdiag(Q_max_i,Q_max_i,Q_max_i);

%%%%%% Decide parameter %%%%%
obs_time = 2; % observation time /second
k_step = obs_time/delta; % time step
r_sensor = 7; %  # of sensors to be chosen

mat_L = gen_mat_L(k_step);

load P_CEKF.mat;
pi_init = P_CEKF{floor(51/delta)+1};
% pi_init = P_CEKF{floor(5/delta)+1};
QCell = repmat({Q_bar},1,k_step);
cov_Z = blkdiag(pi_init,QCell{:});

% R_mea = blkdiag(0.01,0.01);
R_mea = blkdiag(0.01*1e4,0.01*1e4);
VCell = repmat({R_mea},1,num_sensor*(k_step+1));
cov_V = blkdiag(VCell{:});

CCell = repmat({H_col},1,k_step+1);
C_0k = blkdiag(CCell{:});

f_mat_O = @(Sensors)gen_mat_O(Sensors,H_col,k_step);
f_cov_z = @(S) cov_Z - cov_Z*f_mat_O(S)'*(f_mat_O(S)*cov_Z*f_mat_O(S)' + cov_V)^-1*f_mat_O(S)*cov_Z;
f_H = @(S) logdet(f_cov_z(S));
% f_H = @(S) trace(f_cov_z(S));

num_remain = r_sensor;
set_selected = [];
set_possible = 1:num_sensor;
while(num_remain > 0)
    cur_gain = -1e6;
    cur_best_sensor = -1;
    for i = 1:num_sensor
        if isempty(intersect(set_selected,i)) % if sensor i is not selected yet
            tmp = [set_selected i];
            tmp_gain = f_H(set_selected) - f_H(tmp);
            if tmp_gain > cur_gain % if a better cost is achieved
                cur_gain = tmp_gain;
                cur_best_sensor = i;
            end
        end
    end
    if cur_best_sensor ~= -1
        set_selected = [set_selected cur_best_sensor];
    end
    num_remain = num_remain - 1;
end
set_selected
% f_H(cur_sensor_add)

%%%%% customized function %%%%%
function mat_L = gen_mat_L(k)
mat_L = [];
for j = 1:k+1 % how many eye(6) in this row partition
    cur_row = zeros(6,6*(k+1));
    for i = 1:j
        cur_row(1:6,6*i-5:6*i) = eye(6);
    end
    mat_L = [mat_L;cur_row];
end
end

function output = logdet(X)
output = log(det(X));
if isinf(output)
    output = -1e5;
end
end

function mat_O = gen_mat_O(on_sensors,mat_C,k)
len_m = length(mat_C);  % dimension of measurement vector
mat_S = zeros(len_m);   % initialize by zeros matrix

for i=1:length(on_sensors) % modify the selection matrix mat_S by list of turned-on sensors
    pos = on_sensors(i);
    mat_S(2*pos-1:2*pos,2*pos-1:2*pos) = eye(2);
end
mat_SC = mat_S*mat_C;
mat_SC = repmat({mat_SC},1,k+1);
mat_SC = blkdiag(mat_SC{:});
mat_L = gen_mat_L(k);

mat_O = mat_SC*mat_L;

end