%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%created by Solmaz S. Kia%%%%%%%%%%%%%%%
%%%%%%%%%%Last Revised July 2015%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%UCI%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%modified by Qi Yan%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Last Revised August 2018%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
% clear;
% clear classes; % clear all persistent variables
clear all
close all;
addpath('initialization');
addpath('iteration');
addpath('visualization');
rng(1,'twister');
NoiseGen(); % generate noise once before running simulation
%%%%%%%%%%%%%%%%%%%%%%k%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
%missed_onoff = 0: Every one updates
%missed_onoff = 1: Only two robots making relative measurement update
%missed_onoff = 2: No CL
Missed_on_off = 0;

cross_on = 1;

color_plot = 'r';
if Missed_on_off == 0 && cross_on==1
    color_plot = 'r';
end
if Missed_on_off == 0 && cross_on==0
  color_plot = 'b';
end  

if Missed_on_off==2 
  color_plot = 'k';
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%Iteration Settings%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[delta,k_f,~,~,~,RelMea_Table,UpdateOrder,OPT_para] = IterationInit(Missed_on_off);
% delta                 : stepsize
% k_f                   : total number of steps
% RelMea_Table          : relative measurements occurrence table
% UpdateOrder           : update order for sequential updating
% OPT_para              : parameters related to optimization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%Robot Prameters%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X0,X_hat0,P0,N] = RobotInit();
% X0     : initial true poses
% X_hat0 : initial estimated poses and covariances
% P0     : initial collective covariances
% N      : number of robots

n_x=3*N; %%total number of states
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = cell(k_f,N);          % true states
X_hat_DR = cell(k_f,N);   % states by dead reckoning
X_hat_CEKF = cell(k_f,N); % states by standard EKF
X_hat_TEKF = cell(k_f,N); % staets by threshold EKF
X_hat_SATEKF = cell(k_f,N); % states by server-assisted TEKF
X_hat_CUKF = cell(k_f,N); % states by CL method with CD
P_DR = cell(k_f,1);           % collective covariance for dead reckoning
P_CEKF = cell(k_f,1);         % collective covariance matrix for standard EKF
P_TEKF = cell(k_f,1);         % collective covariance matrix for threshold EKF
P_CUKF = cell(k_f,1);         % collective covariance matrix with CD
P_single_SATEKF = cell(k_f,N); % robot-wise covariance matrix for SA TEKF
Phi_SATEKF = cell(k_f,N);     % robot-wise Phi covariance for SA TEKF
% P_single_CD = cell(k_f,N);  % robot-wise covariacne matrix for CD
flag_robot = cell(N,1); % flags for each robot observation
mea_schedule = cell(k_f,1);

% load initial information when k = 0
for i=1:N
    X{1,i} = X0{i};
    X_hat_DR{1,i} = X_hat0{i};
    X_hat_CEKF{1,i} = X_hat0{i};
    X_hat_TEKF{1,i} = X_hat0{i};
    X_hat_SATEKF{1,i} = X_hat0{i};
    X_hat_CUKF{1,i} = X_hat0{i};
    flag_robot{i} = 1;
end
P_DR{1} = P0;  P_CEKF{1} = P0; P_TEKF{1} = P0;
P_CUKF{1} = P0;
for i=1:N
    P_single_SATEKF{1,i} = P0(2*i-1:2*i,2*i-1:2*i);
    Phi_SATEKF{1,i} = eye(2);
end
% clear X0 X_hat0 P0;

load noise_profile.mat;
flag_mea = 0; % flag for measurement (relative or absolute)
flag_seq = 0; % flag for sequantial updating
j_mea = 1;    % times of measurement period
count_com_TEKF = zeros(k_f,2); count_com_SATEKF = zeros(k_f,2);
Three_sigma_DR = zeros(2,k_f,N);
Three_sigma_CEKF = zeros(2,k_f,N); Three_sigma_TEKF = zeros(2,k_f,N);
Three_sigma_SATEKF = zeros(2,k_f,N); Three_sigma_CUKF = zeros(2,k_f,N); % collective 3 sigma boundaries
phi_true = zeros(k_f,N); phi_mea = zeros(k_f,N);
phi_true(1,:) = rand(1,N)*5/57.3;
% phi_true(1,:) = [0 pi/2 pi pi/3 pi/6];
phi_mea = phi_true; % initial orientations
index_noise_extero = 1;
% Three_sigma_robotx_ICI_EKF = zeros(k_f,N); Three_sigma_roboty_ICI_EKF = zeros(k_f,N);
%% Iteration
for k=1:k_f
    disp([num2str(k),' out of ',num2str(k_f)]);
    for j=1:N
        Three_sigma_DR(1,k,j) = 3*sqrt(P_DR{k}(2*j-1,2*j-1));
        Three_sigma_DR(2,k,j) = 3*sqrt(P_DR{k}(2*j,2*j));
        Three_sigma_CEKF(1,k,j) = 3*sqrt(P_CEKF{k}(2*j-1,2*j-1));
        Three_sigma_CEKF(2,k,j) = 3*sqrt(P_CEKF{k}(2*j,2*j));
        Three_sigma_TEKF(1,k,j) = 3*sqrt(P_TEKF{k}(2*j-1,2*j-1));
        Three_sigma_TEKF(2,k,j) = 3*sqrt(P_TEKF{k}(2*j,2*j));
        Three_sigma_SATEKF(1,k,j) = 3*sqrt(P_single_SATEKF{k,j}(1,1));
        Three_sigma_SATEKF(2,k,j) = 3*sqrt(P_single_SATEKF{k,j}(2,2));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % True robot model (pose propagation with noise)
    [X(k+1,:),phi_true(k+1,:),phi_mea(k+1,:),sigma_V,sigma_phi] = TrueRobotModel(X(k,:),phi_true(k,:),phi_mea(k,:),k,Noise_propa_v(k,:),Noise_propa_w(k,:));
    
    % dead reckoning propagation for pose and covariance
    [X_hatProp_DR,P_Prop_DR] = Prop_DR(X_hat_DR(k,:),phi_mea(k+1,:),P_DR{k},k,sigma_V,sigma_phi);
    [X_hatProp_CEKF,P_Prop_CEKF] = Prop_CEKF(X_hat_CEKF(k,:),phi_mea(k+1,:),P_CEKF{k},k,sigma_V,sigma_phi);
    [X_hatProp_TEKF,P_Prop_TEKF] = Prop_CEKF(X_hat_TEKF(k,:),phi_mea(k+1,:),P_TEKF{k},k,sigma_V,sigma_phi);
    [X_hatProp_SATEKF,P_Prop_SATEKF,Phi_Prop_SATEKF,G_Prop_SATEKF] = Prop_SAEKF(X_hat_SATEKF(k,:),phi_mea(k+1,:),P_single_SATEKF(k,:),Phi_SATEKF(k,:),k,sigma_V,sigma_phi);
    X_hat_DR(k+1,:) = X_hatProp_DR; P_DR{k+1} = P_Prop_DR;
    X_hat_CEKF(k+1,:) = X_hatProp_CEKF; P_CEKF{k+1} = P_Prop_CEKF;
    X_hat_TEKF(k+1,:) = X_hatProp_TEKF; P_TEKF{k+1} = P_Prop_TEKF;
    X_hat_SATEKF(k+1,:) = X_hatProp_SATEKF; P_single_SATEKF(k+1,:) = P_Prop_SATEKF;
    Phi_SATEKF(k+1,:) = Phi_Prop_SATEKF;
    count_com_TEKF(k+1,:) = count_com_TEKF(k,:);
    count_com_SATEKF(k+1,:) = count_com_SATEKF(k,:);

    % declaring the master and landmark at the beginning of measurement period
    if RelMea_Table{1,j_mea}(1,1) == k % if it is the starting moment of a measurement period
        cur_start_k = RelMea_Table{1,j_mea}(1,1);
        a_set = RelMea_Table{1,j_mea}(2:end,1); % master robot identities
        b_set = RelMea_Table{1,j_mea}(2:end,2); % landmark robot identities
        size_RelMea_Table = size(RelMea_Table{1,j_mea});

        % elements in MeaMat: 
        % update sequence - robot_a - robot_b - robot_missed
        MeaMat = zeros(length(a_set),size_RelMea_Table(2)+1);
        for ii=1:length(a_set)
            MeaMat(ii,1) = find(UpdateOrder == a_set(ii));
        end
        MeaMat(:,2:end) = RelMea_Table{1,j_mea}(2:end,1:end);
        % Sorting the order of the update order based on the prespecified guidline saved in UpdateOrder
        MeaMat = sortrows(MeaMat,1);
        
        % Create constant measurement matrix
        cur_H_o = cell(N,1);
        cur_R_bar = cell(N,1);
        for ii_a = 1:length(a_set)
            ii_b = b_set(ii_a);
            cur_H_o{ii_a} = zeros(2,2*N);
            cur_H_o{ii_a}(:,2*a_set(ii_a)-1:2*a_set(ii_a)) = -eye(2);
            cur_H_o{ii_a}(:,2*b_set(ii_a)-1:2*b_set(ii_a)) = eye(2);
            [R_rel,R_abs] = ExteroVar();
            if ii_a ~= ii_b
                cur_R_bar{ii_a} = R_rel{a_set(ii_a)}(1:2,1:2);
            elseif ii_a == ii_b
                cur_R_bar{ii_a} = R_abs{a_set(ii_a)}(1:2,1:2);
            end
        end
        cur_H_o = cell2mat(cur_H_o);
        cur_R_bar = blkdiag(cur_R_bar{:});
        
        %%%%%%%%% measurement-related flag settings %%%%%%%%%
%         if size_RelMea_Table(1) > 2 % see if sequential measurement is on
%             flag_seq = 1;
%         end
        flag_seq = 1;
        flag_mea=1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%if there is a reletivemeasurement%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flag_mea == 1
        if mod(k-cur_start_k,OPT_para(1)) == 0
            flag_update = 1;
        else
            flag_update = 0;
        end
        
        %%%%% optimize for measurement scheduling by changing MeaMat %%%%%
        if mod(k-cur_start_k,0.5/delta) == 0 % change topology per some seconds
            P_col = blkdiag(P_single_SATEKF{k,:});
            cov_Z = blkdiag(P_col,G_Prop_SATEKF);
            cell_C = cell(1,N); [R_rel,~] = ExteroVar();
            mat_rotation = @(phi) [cos(phi) -sin(phi);sin(phi) cos(phi)]';
            ab_set = zeros(10,2);
            ii_sensor = 1;
            for i=1:N  % for each robot
                cell_C{i} = zeros(2*(N-1),2*N);
                for j=1:N-1
                    cell_C{i}(2*j-1:2*j,2*i-1:2*i) = mat_rotation(phi_mea(k,i))*-eye(2);
                end
                cur_row=1;
                for j=[1:i-1,i+1:N]
                    cell_C{i}(2*cur_row-1:2*cur_row,2*j-1:2*j) = mat_rotation(phi_mea(k,i))*eye(2);
                    cur_row=cur_row+1;
                    ab_set(ii_sensor,:) = [i,j];
                    ii_sensor = ii_sensor+1;
                end
            end
            mat_C = vertcat(cell_C{:});
            fun_mat_O = @(S) OPT_get_O(S,mat_C);
            cov_V = repmat(R_rel,1,2*(N-1));
            cov_V = blkdiag(cov_V{:});
            cov_Z_post = @(S)cov_Z - cov_Z*fun_mat_O(S)'*(fun_mat_O(S)*cov_Z*fun_mat_O(S)' + cov_V)^-1*fun_mat_O(S)*cov_Z;
            f_logdet = @(X)log(det(X));
            f_H = @(S)f_logdet(cov_Z_post(S));
            
            num_sensor = N*(N-1);
            r_sensor = N-2;
            num_remain = r_sensor;
%             set_random = randi(num_sensor,1,r_sensor);
%             set_random = [1 6 11 16 17];
            set_random = 1:num_sensor;
            f_H_random = f_H(set_random);
            set_selected = [];
            set_possible = 1:num_sensor;
            cur_H=1e6;
            while(num_remain > 0 && cur_H > f_H_random)
%             while(num_remain > 0)
                cur_gain = -1e6;
                cur_best_sensor = -1;
                for i = 1:num_sensor
                    if isempty(intersect(set_selected,i)) % if sensor i is not selected yet
                        tmp = [set_selected i]; % select sensor i
                        tmp_gain = f_H(set_selected) - f_H(tmp); % calculate the gain of selecting sensor i
                        if tmp_gain > cur_gain % if a better gain is achieved
                            cur_gain = tmp_gain;
                            cur_best_sensor = i;
                        end
                    end
                end
                if cur_best_sensor ~= -1
                    set_selected = [set_selected cur_best_sensor];
                end
                cur_H = f_H(set_selected);
                num_remain = num_remain - 1;
            end
                    
            MeaMat = []; MeaMat_random = [];
            for i=1:length(set_selected)
                MeaMat = [MeaMat; i ab_set(set_selected(i),:)];
            end
            for i=1:length(set_random)
                MeaMat_random = [MeaMat_random;i ab_set(set_random(i),:)];
            end
            a_set = MeaMat(:,2);
            b_set = MeaMat(:,3);
            a_set_random = MeaMat_random(:,2);
            b_set_random = MeaMat_random(:,3);
        end
        %%%%% end %%%%%
        %%%%% optimize for measurement scheduling by changing MeaMat %%%%%
        
        % after opt decision, decide whether to measurement or not for each
        % robot
        if k > cur_start_k
            if mod(k-cur_start_k,X_opt(1)) == 0
                flag_robot{1} = 1;
            else
                flag_robot{1} = 0;
            end
            if mod(k-cur_start_k,X_opt(2)) == 0
                flag_robot{2} = 1;
            else
                flag_robot{2} = 0;
            end
            if mod(k-cur_start_k,X_opt(3)) == 0
                flag_robot{3} = 1;
            else
                flag_robot{3} = 0;
            end
            if mod(k-cur_start_k,update_interval) == 0
                flag_update = 1; % update at server-side
            else
                flag_update = 0;
            end
        end
        
        % decide the optimal observation rate and communication rate
        % *after the starting time step* %
        if k == cur_start_k
%             [V,~] = TrajGen(k);
%             [sigma_V,~] = ProprioVar(V,k);
%             Q_max = zeros(N,1); Q_bar = cell(N,1);
%             for ii = 1:length(V)
%                 Q_max(ii) = sigma_V(ii)^2*V(ii);
%                 Q_bar{ii} = Q_max(ii)*eye(2);
%             end
%             Q_bar = blkdiag(Q_bar{:});
%             
% %             Pi_ss = @(W)dare(eye(2*N),cur_H_o',delta^2*Q_bar,W,zeros(2*N),eye(2*N));
%             Pi_ss = @(W)dare(eye(2*N),cur_H_o',delta^2*Q_bar,W);
% %             Pi_ss = @(W)care(zeros(2*N),chol(cur_H_o'*W^-1*cur_H_o,'lower'),delta^2*Q_bar);
%             mat_f = @(f) blkdiag(1/delta/f(1)*eye(2),1/delta/f(2)*eye(2),1/delta/f(3)*eye(2));
%             Pi_bar = @(f,Z)Pi_ss(cur_R_bar*mat_f(f) + Z^-1);
%             f_rate = @(f,Z) 1 - 1/sqrt(det(eye(2*N)+(cur_H_o*Pi_bar(f,Z)*cur_H_o' + cur_R_bar*mat_f(f))*Z));
%             R_star = @(f,Z) (f_rate(f,Z)*(cur_R_bar*mat_f(f))^-1 + (1-f_rate(f,Z))*(cur_R_bar*mat_f(f) + Z^-1)^-1)^-1;
%             Pi_exp = @(f,Z) Pi_ss(R_star(f,Z));
%             mu_o = 40; mu_c = 100;
%             f_cost = @(X) mu_o*sum(X(1:3)) + mu_c*f_rate(X(1:3),X(4)*eye(2*N))*max(X(1:3));
%             
%             f_cost_int = @(X) f_cost([1/delta*[1/X(1) 1/X(2) 1/X(3)] X(4)]);
%             f_rate_int = @(X) f_rate(1/delta*[1/X(1) 1/X(2) 1/X(3)],X(4)*eye(2*N));
%             X_init = OPT_para;
%             Pi_exp_init = Pi_exp(1/delta*X_init(1:3).^-1,X_init(4)*eye(2*N));
%             Pi_bar_init = Pi_bar(1/delta*X_init(1:3).^-1,X_init(4)*eye(2*N));
%             cost_init = f_cost_int(X_init);
%             %                 f_con = @(X) trace(Pi_exp(X(1:3).^-1,X(4)*eye(2*N))) - trace(Pi_exp_init);
%             f_con = @(X) trace(Pi_bar(1/delta*X(1:3).^-1,X(4)*eye(2*N))) - 0.95*trace(Pi_bar_init);
%             f_nlcon = @(X)deal(f_con(X),[]);
%             [X_opt,cost_opt] = ga(f_cost_int,4,[],[],[],[],[1 1 1 1],[10 10 10 100],f_nlcon,[1 2 3 4]);
%             update_interval = min(X_opt(1:3));
%             coef_Z = X_opt(end);
%             Pi_bar_opt = Pi_bar(1/delta*X_opt(1:3).^-1,X_opt(4)*eye(2*N));
%             ratio_cost = cost_opt/cost_init
%             ratio_tr = trace(Pi_bar_opt)/trace(Pi_bar_init)
%             clear V sigma_V;
            X_opt = [1 1 1 1e9];
            update_interval = min(X_opt(1:3));
            coef_Z = X_opt(end);
        end
        
        mea_schedule{k} = MeaMat;
        %ii_a - relative measurment order
        for ii_a = 1:length(a_set)
            cur_a = MeaMat(ii_a,2);
            cur_b = MeaMat(ii_a,3);
            if size_RelMea_Table(2) > 2 % more than 2 columns
                miss_set = MeaMat(ii_a,4:end);
            else
                miss_set = 0;
            end
            % double check the measurement period
            if (k>=RelMea_Table{1,j_mea}(1,1) && k<=RelMea_Table{1,j_mea}(1,2))
                if flag_seq == 1 && ii_a > 1 % not the first measurement at time step k
                    [K_CEKF,S_ab_CEKF,r_a_CEKF] = Update_CEKF(cur_a,cur_b,X{k+1,cur_a},X{k+1,cur_b},phi_mea(k+1,:),X_hat_CEKF{k+1,cur_a},X_hat_CEKF{k+1,cur_b},P_CEKF{k+1},k,Noise_extero(index_noise_extero,:));
                    if flag_robot{cur_a} == 1
                        [r_a_bar,Gamma_col,DP_SATEKF,param_SATEKF] = Update_SATEKF(cur_a,cur_b,X{k+1,cur_a},X{k+1,cur_b},phi_mea(k+1,:),X_hat_SATEKF{k+1,cur_a},P_single_SATEKF{k+1,cur_a},Phi_SATEKF{k+1,cur_a},X_hat_SATEKF{k+1,cur_b},P_single_SATEKF{k+1,cur_b},Phi_SATEKF{k+1,cur_b},k,Noise_extero(index_noise_extero,:),coef_Z);
                        count_com_SATEKF(k+1,:) = count_com_SATEKF(k+1,:) + 1;
                    end
                    index_noise_extero = index_noise_extero + 1;
                elseif ii_a == 1 % the first measurement at time step k
                    [K_CEKF,S_ab_CEKF,r_a_CEKF] = Update_CEKF(cur_a,cur_b,X{k+1,cur_a},X{k+1,cur_b},phi_mea(k+1,:),X_hat_CEKF{k+1,cur_a},X_hat_CEKF{k+1,cur_b},P_CEKF{k+1},k,Noise_extero(index_noise_extero,:));
                    if flag_robot{cur_a} == 1
                        [r_a_bar,Gamma_col,DP_SATEKF,param_SATEKF] = Update_SATEKF(cur_a,cur_b,X{k+1,cur_a},X{k+1,cur_b},phi_mea(k+1,:),X_hat_SATEKF{k+1,cur_a},P_single_SATEKF{k+1,cur_a},Phi_SATEKF{k+1,cur_a},X_hat_SATEKF{k+1,cur_b},P_single_SATEKF{k+1,cur_b},Phi_SATEKF{k+1,cur_b},k,Noise_extero(index_noise_extero,:),coef_Z);
                        count_com_SATEKF(k+1,:) = count_com_SATEKF(k+1,:) + 1;
                    end
                    index_noise_extero = index_noise_extero + 1;
                end
            end
            % update for SATEKF
            if flag_robot{cur_a} == 1 % robot-to-server update is allowed
                for i = 1:N
                    if sum(i == miss_set) == 0 % robot i doesn't miss update
                        X_hat_SATEKF{k+1,i} = X_hat_SATEKF{k+1,i} + Phi_SATEKF{k+1,i}*Gamma_col{i}*r_a_bar;
                        P_single_SATEKF{k+1,i} = P_single_SATEKF{k+1,i} - Phi_SATEKF{k+1,i}*(Gamma_col{i}*Gamma_col{i}'+DP_SATEKF{i})*Phi_SATEKF{k+1,i}';
                    end
                end
                if param_SATEKF(1) == 0
                    count_com_SATEKF(k+1,1) = count_com_SATEKF(k+1,1) - 1;
                end
            end
            
            % update for benchmark CEKF
            KK_CEKF = cell(1,N);
            for i=1:N
                KK_CEKF{i} = K_CEKF(2*i-1:2*i,:);
                if sum(i == miss_set) == 0 % robot i doesn't miss update
                    X_hat_CEKF{k+1,i} = X_hat_CEKF{k+1,i} + KK_CEKF{i}*r_a_CEKF;
                end
            end
            P_CEKF{k+1} = P_CEKF{k+1} - K_CEKF*S_ab_CEKF*K_CEKF';
        end % end of ii_a with respect to optimized MeaMat
        
        %ii_a - relative measurment order
        for ii_a = 1:length(a_set_random)
            cur_a = MeaMat_random(ii_a,2);
            cur_b = MeaMat_random(ii_a,3);
            if size_RelMea_Table(2) > 2 % more than 2 columns
                miss_set = MeaMat_random(ii_a,4:end);
            else
                miss_set = 0;
            end
            % double check the measurement period
            if (k>=RelMea_Table{1,j_mea}(1,1) && k<=RelMea_Table{1,j_mea}(1,2))
                if flag_seq == 1 && ii_a > 1 % not the first measurement at time step k
                    if flag_robot{cur_a} == 1
                        [K_TEKF,S_ab_TEKF,r_a_TEKF,param_TEKF] = Update_TEKF(cur_a,cur_b,X{k+1,cur_a},X{k+1,cur_b},phi_mea(k+1,:),X_hat_TEKF{k+1,cur_a},X_hat_TEKF{k+1,cur_b},P_TEKF{k+1},k,Noise_extero(index_noise_extero,:),coef_Z);
                        count_com_TEKF(k+1,:) = count_com_TEKF(k+1,:) + 1;
                    end
                    index_noise_extero = index_noise_extero + 1;
                elseif ii_a == 1 % the first measurement at time step k
                    if flag_robot{cur_a} == 1
                        [K_TEKF,S_ab_TEKF,r_a_TEKF,param_TEKF] = Update_TEKF(cur_a,cur_b,X{k+1,cur_a},X{k+1,cur_b},phi_mea(k+1,:),X_hat_TEKF{k+1,cur_a},X_hat_TEKF{k+1,cur_b},P_TEKF{k+1},k,Noise_extero(index_noise_extero,:),coef_Z);
                        count_com_TEKF(k+1,:) = count_com_TEKF(k,:) + 1;
                    end
                    index_noise_extero = index_noise_extero + 1;
                end
            end
            % update for TEKF
            if flag_robot{cur_a} == 1 % robot-to-server update is allowed
                KK_TEKF = cell(1,N);
                for i = 1:N
                    KK_TEKF{i} = K_TEKF(2*i-1:2*i,:);
                    if sum(i == miss_set) == 0 % robot i doesn't miss update
                        X_hat_TEKF{k+1,i} = X_hat_TEKF{k+1,i} + param_TEKF(1)*KK_TEKF{i}*r_a_TEKF;
                    end
                end
                P_TEKF{k+1} = P_TEKF{k+1} - param_TEKF(2)*K_TEKF*S_ab_TEKF*K_TEKF';
                if param_TEKF(1) == 0 % under threshold
                    count_com_TEKF(k+1,1) = count_com_TEKF(k+1,1) - 1;
                end
            end
        end % end of ii_a with respect to randomized MeaMat
        
        % at the end of each measurement period, set flag variables
        if RelMea_Table{1,j_mea}(1,2) == k
            j_mea = j_mea+1;
            flag_mea = 0;
            clear MeaMat a_set b_set miss_set;
        end
    end % end of relative flag
    
end  % end of k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve DARE for upper bound
% H_i = zeros(2,N*2);
% H_i(:,1:2) = eye(2); H_i(:,3:4) = -eye(2);
% care_B = H_i'/sqrt(delta);
% Q_max = [0.05*0.3 0; 0 0.05^2*0.3];
% Q_max = blkdiag(Q_max,Q_max,Q_max,Q_max,Q_max);
% care_Q = delta*Q_max;
% care_R = blkdiag(0.04^2,0.04^2);
% [X_dare,~,~] = dare(zeros(10),care_B,care_Q,care_R,zeros(10,2),eye(10));
%% Visualization
XX = zeros(2,k_f+1,N); XX_hat_DR = zeros(2,k_f+1,N); 
XX_hat_CEKF = zeros(2,k_f+1,N); XX_hat_TEKF = zeros(2,k_f+1,N);
XX_hat_SATEKF = zeros(2,k_f+1,N); XX_hat_CUKF = zeros(2,k_f+1,N);
for i=1:N
    XX(:,:,i)=cell2mat(X(:,i)');
    XX_hat_DR(:,:,i) = cell2mat(X_hat_DR(:,i)');
    XX_hat_CEKF(:,:,i) = cell2mat(X_hat_CEKF(:,i)');
    XX_hat_TEKF(:,:,i) = cell2mat(X_hat_TEKF(:,i)');
    XX_hat_SATEKF(:,:,i) = cell2mat(X_hat_SATEKF(:,i)');
end
k=k_f+1;
for j=1:N
    Three_sigma_DR(1,k,j) = 3*sqrt(P_DR{k}(2*j-1,2*j-1));
    Three_sigma_DR(2,k,j) = 3*sqrt(P_DR{k}(2*j,2*j));
    Three_sigma_CEKF(1,k,j) = 3*sqrt(P_CEKF{k}(2*j-1,2*j-1));
    Three_sigma_CEKF(2,k,j) = 3*sqrt(P_CEKF{k}(2*j,2*j));
    Three_sigma_TEKF(1,k,j) = 3*sqrt(P_TEKF{k}(2*j-1,2*j-1));
    Three_sigma_TEKF(2,k,j) = 3*sqrt(P_TEKF{k}(2*j,2*j));
    Three_sigma_SATEKF(1,k,j) = 3*sqrt(P_single_SATEKF{k,j}(1,1));
    Three_sigma_SATEKF(2,k,j) = 3*sqrt(P_single_SATEKF{k,j}(2,2));
end
NEES_DR = zeros(1,k_f+1); NEES_CEKF = zeros(1,k_f+1);
NEES_TEKF = zeros(1,k_f+1);  NEES_SATEKF = zeros(1,k_f+1);
NEES_CUKF = zeros(1,k_f+1);
Tr_DR = NEES_DR; Tr_CEKF = NEES_CEKF;
Tr_TEKF = NEES_CEKF; Tr_SATEKF = NEES_SATEKF;
Tr_CUKF = NEES_CUKF;
RMSE_DR = zeros(k_f+1,N); RMSE_CEKF = zeros(k_f+1,N);
RMSE_TEKF = zeros(k_f+1,N); RMSE_SATEKF = zeros(k_f+1,N);
RMSE_CUKF = zeros(k_f+1,N);
for j=1:length(NEES_CEKF)
    X_wave = cat(1,X{j,:}) - cat(1,X_hat_DR{j,:});
    NEES_DR(j) = X_wave'*P_DR{j}^-1*X_wave/n_x;
    Tr_DR(j) = trace(P_DR{j});
    tmp = sqrt(reshape(X_wave.^2,2,N));
    RMSE_DR(j,:) = sum(tmp(1:2,:),1);

    X_wave = cat(1,X{j,:}) - cat(1,X_hat_CEKF{j,:});
    NEES_CEKF(j) = X_wave'*P_CEKF{j}^-1*X_wave/n_x;
    Tr_CEKF(j) = trace(P_CEKF{j});
    tmp = sqrt(reshape(X_wave.^2,2,N));
    RMSE_CEKF(j,:) = sum(tmp(1:2,:),1);
    
    X_wave = cat(1,X{j,:}) - cat(1,X_hat_TEKF{j,:});
    NEES_TEKF(j) = X_wave'*P_TEKF{j}^-1*X_wave/n_x;
    Tr_TEKF(j) = trace(P_TEKF{j});
    tmp = sqrt(reshape(X_wave.^2,2,N));
    RMSE_TEKF(j,:) = sum(tmp(1:2,:),1);
    
    X_wave = cat(1,X{j,:}) - cat(1,X_hat_SATEKF{j,:});
    NEES_SATEKF(j) = X_wave'*blkdiag(P_single_SATEKF{j,:})^-1*X_wave/n_x;
    Tr_SATEKF(j) = trace(blkdiag(P_single_SATEKF{j,:}));
    tmp = sqrt(reshape(X_wave.^2,2,N));
    RMSE_SATEKF(j,:) = sum(tmp(1:2,:),1);
end
Tk = (0:delta:k_f*delta);
t_f = k_f*delta;

disp(['TEKF Number of communications: 1)with threshold ',num2str(count_com_TEKF(end,1)),' ,2)without threshold ',num2str(count_com_TEKF(end,2))]);
disp(['SATEKF Number of communications: 1)with threshold ',num2str(count_com_SATEKF(end,1)),' ,2)without threshold ',num2str(count_com_SATEKF(end,2))]);
save CL_data.mat t_f delta Tk k_f N RelMea_Table color_plot Missed_on_off
save CL_data.mat XX XX_hat_DR Three_sigma_DR RMSE_DR -append;
save CL_data.mat XX_hat_CEKF Three_sigma_CEKF RMSE_CEKF -append;
save CL_data.mat XX_hat_TEKF Three_sigma_TEKF RMSE_TEKF -append;
% save CL_data.mat XX_hat_CUKF Three_sigma_CUKF RMSE_CUKF -append;
% XX_hat_CD XX_hat_ICI XX_hat_ICI_EKF
% save CL_data.mat Three_sigma_robotx_CD Three_sigma_roboty_CD NEES_CD NEES_ICI NEES_ICI_EKF -append;
% save CL_data.mat Three_sigma_robotx_ICI_EKF Three_sigma_roboty_ICI_EKF -append;
% save CL_data.mat Tr_CEKF Tr_CD Tr_ICI Tr_ICI_EKF -append;
% save CL_data.mat
save P_CEKF.mat P_CEKF;
Plot_Results;
% Figures_Plot; % get figure windows
% Traj_Plot; % plot results
%%
toc;