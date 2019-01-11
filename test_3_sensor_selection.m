clear
close all

% Construct observation matrix
H = cell(7,1);
H{1} = [-eye(2) eye(2) zeros(2)];       % 1->2
H{2} = [-eye(2) zeros(2) eye(2)];       % 1->3
H{3} = [eye(2) -eye(2) zeros(2)];       % 2->1
H{4} = [zeros(2) -eye(2) eye(2)];       % 2->3
H{5} = [zeros(2) eye(2) -eye(2)];       % 3->2
H{6} = [eye(2) zeros(2) -eye(2)];       % 3->1
H{7} = [zeros(2) zeros(2) -eye(2)];     % 3->3
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

base_I = zeros(len_m,1);
f_mat_I = @(i) gen_mat_I(i,len_m,k_step);
f_mat_M = @(i) C_0k'*f_mat_I(i)*cov_V^-1*f_mat_I(i)*C_0k;
f_sensor_add = @(i) mat_L'*f_mat_M(i)*mat_L;
f_H = @(X) -logdet(X + cov_Z^-1); % logdet of posterior cov_Z
% f_H = @(X) (-logdet(X + cov_Z^-1)); % -logdet of posterior cov_Z
f_H_delta = @(X,i) f_H(X+f_sensor_add(i)) - f_H(X); % negative delta
f_H_gain = @(X,i) f_H(X) - f_H(X+f_sensor_add(i)); % positive gain

num_remain = r_sensor;
set_selected = [];
set_possible = 1:num_sensor;
cur_sensor_init = zeros(length(f_sensor_add(1)));
cur_sensor_add = cur_sensor_init;
while(num_remain > 0)
    cur_gain = -1e6;
    cur_best_sensor = -1;
    for i = 0:num_sensor
        if isempty(intersect(set_selected,i)) % if sensor i is not selected yet
            if f_H_gain(cur_sensor_add,i) > cur_gain % if a better cost is achieved
                cur_gain = f_H_gain(cur_sensor_add,i);
                cur_best_sensor = i;
            end
        end
    end
    if cur_best_sensor ~= -1
        cur_sensor_add = cur_sensor_add + f_sensor_add(cur_best_sensor);
    end
    set_selected = [set_selected cur_best_sensor];
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

function mat_I = gen_mat_I(i,n,k)
if 2*i > n || i < 0
    disp('Bad input for i!');
    mat_I = zeros(n);
    return;
end
if i==0
    mat_I = zeros(length((k+1)*n));
    return;
end
base = zeros(n,1);
base(2*i-1) = 1;
base(2*i) = 1;
base = diag(base);
base_cell = repmat({base},1,k+1);
mat_I = blkdiag(base_cell{:});
end

function output = logdet(X)
output = log(det(X));
if isinf(output)
    output = -1e5;
end
end