%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%created by Solmaz S. Kia%%%%%%%%%%%%%%%
%%%%%%%%%%Last Revised July 2015%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%UCI%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%modified by Qi Yan%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Last Revised August 2018%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
clear classes; % clear all persistent variables
% close all;
addpath('initialization');
addpath('iteration');
addpath('visualization');
Monte_Carlo_Profile_Gen(); % generate noise once before running simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clc
% clear X_hat XX_hat P_bar P P_Prop Phi a_set b_set P0;
% clear XX X;
% clear MeaMat Three_sigma_robotx  Three_sigma_roboty;
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
[delta,k_f,~,~,~,RelMea_Table,UpdateOrder] = IterationInit(Missed_on_off);
% delta                 : stepsize
% k_f                   : total number of steps
% RelMea_Table          : relative measurements occurrence table
% UpdateOrder           : update order for sequential updating

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%Robot Prameters%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X0,X_hat0,P0,N] = RobotInit();
% X0     : initial true poses
% X_hat0 : initial estimated poses and covariances
% P0     : initial collective covariances
% N      : number of robots

n_x=3*N; %%total number of states
Nx = 3*ones(1,N);

load velocity_profile_mont.mat
load noise_profile_mont.mat
load AnalysisData.mat

%% Iteration
for Monte_index=1:M
    Monte_index
    W=W_mont{Monte_index};
    Noise_propa_w=Noise_propa_w_mont{Monte_index};
    Noise_propa_v=Noise_propa_v_mont{Monte_index};
    
    Noise_extero=[Noise_extero_x_mont{Monte_index} Noise_extero_y_mont{Monte_index} ...
        Noise_extero_phi_mont{Monte_index}  Noise_extero_ax_mont{Monte_index}  Noise_extero_ay_mont{Monte_index}];
X = cell(k_f,N);         % true states
X_hat_CD = cell(k_f,N);  % states by CL method with CD
X_hat_SD = cell(k_f,N);  % states by CL method with SD
X_hat_IN = cell(k_f,N);  % states by CL method with IN
P_CD = cell(k_f,1);      % collective covariance matrix with CD
P_SD = cell(k_f,1);      % collective covariance matrix with SD
P_IN = cell(k_f,1);      % collective covariance matrix with IN
P_single = cell(k_f,N);  % robot-wise covariacne matrix (for CD only)

% load initial information when k = 0
for i=1:N
    X{1,i}=X0{i};
    X_hat_CD{1,i}=X_hat0{i};
    X_hat_SD{1,i}=X_hat0{i};
    X_hat_IN{1,i}=X_hat0{i};
end
P_CD{1}=P0; P_SD{1}=P0; P_IN{1}=P0;
P_single{1,1}=P_CD{1}(1:Nx(1),1:Nx(1));
for i=2:N
    P_single{1,i}=P_CD{1}(sum(Nx(1:i-1))+1:sum(Nx(1:i)),sum(Nx(1:i-1))+1:sum(Nx(1:i)));
end
% clear X0 X_hat0 P0;

load noise_profile.mat;
flag_mea = 0; % flag for measurement (relative or absolute)
flag_seq = 0; % flag for sequantial updating
j_mea = 1;    % times of measurement
Three_sigma_robotx = zeros(k_f,N); Three_sigma_roboty = zeros(k_f,N); % collective 3 sigma boundaries
%% Iteration
for k=1:k_f
    for j=1:N
        Three_sigma_robotx(k,j)=3*sqrt(P_single{k,j}(1,1));
        Three_sigma_roboty(k,j)=3*sqrt(P_single{k,j}(2,2));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % True robot model (pose propagation with noise)
    [X(k+1,:),sigma_V,sigma_W] = TrueRobotModel(X(k,:),k,Noise_propa_v(k,:),Noise_propa_w(k,:));
    
    % dead reckoning propagation for pose and covariance
    [X_hatProp_CD,P_Prop_CD,Chi_x_minus_CD,E_x_CD] = UKFCL_Propagation_CD(X_hat_CD(k,:),P_CD{k},k,sigma_V,sigma_W);
    [X_hatProp_SD,P_Prop_SD,Chi_x_minus_SD,E_x_SD] = UKFCL_Propagation_SD(X_hat_SD(k,:),P_SD{k},k,sigma_V,sigma_W);
    [X_hatProp_IN,P_Prop_IN,Chi_x_minus_IN,E_x_IN] = UKFCL_Propagation_IN(X_hat_IN(k,:),P_IN{k},k,sigma_V,sigma_W);
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%if there is no reletivemeasurement%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X_hat_CD(k+1,:) = X_hatProp_CD;
    P_CD{k+1} = 0.5*(P_Prop_CD+P_Prop_CD');
    
    X_hat_SD(k+1,:) = X_hatProp_SD;
    P_SD{k+1} = 0.5*(P_Prop_SD+P_Prop_SD');    
    
    X_hat_IN(k+1,:) = X_hatProp_IN;
    P_IN{k+1} = 0.5*(P_Prop_IN+P_Prop_IN');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%if there is a reletivemeasurement%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % declaring the interim master and landmark at the beginning of
    % measurement period
    if RelMea_Table{1,j_mea}(1,1) == k
        a_set = RelMea_Table{1,j_mea}(2:end,1);
        b_set = RelMea_Table{1,j_mea}(2:end,2);
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
        
        if size_RelMea_Table(1) > 2
            flag_seq = 1;
        end
        flag_mea=1;
    end
      
    if flag_mea==1
        %ii_a %relative measurment
        for ii_a = 1:length(a_set)
            a = MeaMat(ii_a,2);
            b = MeaMat(ii_a,3);
            if size_RelMea_Table(1,2) > 2
                RobotMissedUpdate = MeaMat(ii_a,4:end);
            else
                RobotMissedUpdate = 0;
            end
            
            % double check the measurement period
            if (k>=RelMea_Table{1,j_mea}(1,1) && k<=RelMea_Table{1,j_mea}(1,2))
                if flag_seq == 1 && ii_a > 1 % i>2 measurements, update sigma points(Chi_x and E_x)
                    [Chi_x_CD,E_x_est_CD] = UKFCL_SeqUpdate_CD(X_hat_CD(k+1,:),P_CD{k+1});
                    [K_UKF_CD,S_ab_CD,r_a_CD] = UKFCL_InterimMaster(a,b,Chi_x_CD,E_x_est_CD,X{k+1,a},X{k+1,b},k,Noise_extero(k,:));
                    [Chi_x_SD,E_x_est_SD] = UKFCL_SeqUpdate_SD(X_hat_SD(k+1,:),P_SD{k+1});
                    [K_UKF_SD,S_ab_SD,r_a_SD] = UKFCL_InterimMaster(a,b,Chi_x_SD,E_x_est_SD,X{k+1,a},X{k+1,b},k,Noise_extero(k,:));
                    [Chi_x_IN,E_x_est_IN] = UKFCL_SeqUpdate_IN(X_hat_IN(k+1,:),P_IN{k+1});
                    [K_UKF_IN,S_ab_IN,r_a_IN] = UKFCL_InterimMaster(a,b,Chi_x_IN,E_x_est_IN,X{k+1,a},X{k+1,b},k,Noise_extero(k,:));
                else % first measurement
                    [K_UKF_CD,S_ab_CD,r_a_CD] = UKFCL_InterimMaster(a,b,Chi_x_minus_CD,E_x_CD,X{k+1,a},X{k+1,b},k,Noise_extero(k,:));
                    [K_UKF_SD,S_ab_SD,r_a_SD] = UKFCL_InterimMaster(a,b,Chi_x_minus_SD,E_x_SD,X{k+1,a},X{k+1,b},k,Noise_extero(k,:));
                    [K_UKF_IN,S_ab_IN,r_a_IN] = UKFCL_InterimMaster(a,b,Chi_x_minus_IN,E_x_IN,X{k+1,a},X{k+1,b},k,Noise_extero(k,:));
                end
            end
            
            KK_CD = cell(1,N); KK_SD = cell(1,N); KK_IN = cell(1,N);
            KK_CD{1}=K_UKF_CD(1:Nx(1),:); KK_SD{1}=K_UKF_SD(1:Nx(1),:);
            KK_IN{1}=K_UKF_IN(1:Nx(1),:);
            for i=2:N
                KK_CD{i}=K_UKF_CD(sum(Nx(1:i-1))+1:sum(Nx(1:i)),:);
                KK_SD{i}=K_UKF_SD(sum(Nx(1:i-1))+1:sum(Nx(1:i)),:);
                KK_IN{i}=K_UKF_IN(sum(Nx(1:i-1))+1:sum(Nx(1:i)),:);
            end
            
            % update for states and covariance
            for i=1:N
                if sum(i == RobotMissedUpdate) == 0 % robot i doesn't miss update
                    X_hat_CD{k+1,i}=X_hat_CD{k+1,i}+KK_CD{i}*r_a_CD;
                    X_hat_SD{k+1,i}=X_hat_SD{k+1,i}+KK_SD{i}*r_a_SD;
                    X_hat_IN{k+1,i}=X_hat_IN{k+1,i}+KK_IN{i}*r_a_IN;
                end
            end
            P_CD{k+1}=P_CD{k+1} - 0.5*(K_UKF_CD*S_ab_CD*K_UKF_CD'+(K_UKF_CD*S_ab_CD*K_UKF_CD')');
            P_SD{k+1}=P_SD{k+1} - 0.5*(K_UKF_SD*S_ab_SD*K_UKF_SD'+(K_UKF_SD*S_ab_SD*K_UKF_SD')');
            P_IN{k+1}=P_IN{k+1} - 0.5*(K_UKF_IN*S_ab_IN*K_UKF_IN'+(K_UKF_IN*S_ab_IN*K_UKF_IN')');
        end % end of ii_a
        
        % at the end of each measurement period
        if RelMea_Table{1,j_mea}(1,2) == k
            j_mea = j_mea+1;
            flag_mea=0;
            clear MeaMat a_set b_set;
        end
    end % end of relative flag
    
    P_single{k+1,1} = P_CD{k+1}(1:Nx(1),1:Nx(1));
    for i=2:N
        P_single{k+1,i} = P_CD{k+1}(sum(Nx(1:i-1))+1:sum(Nx(1:i)),sum(Nx(1:i-1))+1:sum(Nx(1:i)));
    end
end  % end of k
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    X_mont{Monte_index} = X;
    X_hat_CD_mont{Monte_index} = X_hat_CD;
    P_CD_mont{Monte_index}= P_CD;
    X_hat_SD_mont{Monte_index} = X_hat_SD;
    P_SD_mont{Monte_index}= P_SD;
    X_hat_IN_mont{Monte_index} = X_hat_IN;
    P_IN_mont{Monte_index}= P_IN;

%     clear X_hat XX_hat P P_Prop  a_set b_set e;
%     clear XX X Noise_propa_w ;
%     clear Noise_propa_v Noise_extero;
%     clear MeasurementMatrix Three_sigma_robotx  Three_sigma_roboty;
end

%% Computation for visualization
% e_pos=zeros(k_f,N);
% RMS_pos=zeros(k_f,N);
% %e_phi=zeros(k_f,N);
% %RMS_phi=zeros(k_f,N);
% for time_step=1:k_f
%     for Monte_index=1:M
%         for UID=1:N
%             error=(X_mont{Monte_index}{time_step,UID}-X_hat_CD_mont{Monte_index}{time_step,UID});
%             error_pos=error(1:2,1);
%             %             error_phi=error(3,1);
%             %             e_pos(time_step,UID)=e_pos(time_step,UID)+error_pos'*inv(P_mont{Monte_index}{time_step}{UID,UID}(1:2,1:2))*error_pos/M;
%             RMS_pos(time_step,UID)=RMS_pos(time_step,UID)+error_pos'*error_pos/M;
%             % e_phi(time_step,UID)=e_phi(time_step,UID)+error_phi'*inv(P_mont{Monte_index}{time_step}{UID,UID}(3,3))*error_phi/M;
%             %RMS_phi(time_step,UID)=RMS_phi(time_step,UID)+error_phi'*error_phi/M;
%             clear error_pos error_phi;
%         end
%     end
% end
% RMS_pos = sqrt(RMS_pos);

disp('Calculating NESS......');
NEES_CD = zeros(1,length(X)); NEES_SD = NEES_CD;
NEES_IN = NEES_CD;
for j=1:length(NEES_CD)
    X_wave = cat(1,X{j,:}) - cat(1,X_hat_CD{j,:});
    NEES_CD(j) = X_wave'*P_CD{j,1}^-1*X_wave/n_x;
    X_wave = cat(1,X{j,:}) - cat(1,X_hat_SD{j,:});
    NEES_SD(j) = X_wave'*P_SD{j,1}^-1*X_wave/n_x;
     X_wave = cat(1,X{j,:}) - cat(1,X_hat_IN{j,:});
    NEES_IN(j) = X_wave'*P_IN{j,1}^-1*X_wave/n_x;
end
Tk = 0:delta:k_f*delta;
save CL_data_mont;

%%
clear
load CL_data_mont
% figure;
% % plot(Tk(1:end-1),e_pos')
% hold on
% % plot(Tk(1:end-1),RMS_pos',color_plot);
% plot(Tk(1:end-1),RMS_pos');
% str_lgd = cell(1,N);
% for i=1:N
%     str_lgd{i} = ['Robot ',num2str(i)];
% end
% text(50,0.5*max(ylim),['M = ',num2str(M)]);
% legend(str_lgd);

% NEES result
figure;
hold on;
plot(Tk,NEES_CD,'linewidth',1);
plot(Tk,NEES_SD,'linewidth',1);
plot(Tk,NEES_IN,'linewidth',1);
line([0 max(xlim)],[1 1],'linewidth',1,'color','red','linestyle','--');
legend('Cholesky Decomposition','Schur Decomposition','IN Iteration','NEES = 1');
title(['Average NEES with M = ',num2str(M)]);
%%
toc;