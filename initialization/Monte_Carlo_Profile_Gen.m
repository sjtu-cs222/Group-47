%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%created by Solmaz S. Kia%%%%%%%%%%%%%%%
%%%%%%%%%%Last Revised April 2014%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%UCSD%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%modified by Qi Yan%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Last Revised August 2018%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Monte_Carlo_Profile_Gen()
%%
%%%%%%%%%%%%%%%%%%%%%%Monte Carlo Run Setup%%%%%%%%%%%%%%%%%%%%%%%%%
% M=50;            %%% # Monte Carlo run
M = 50;            %%% # Monte Carlo run

[~,k_f,~,~,~,~,~] = IterationInit();
length_noise = k_f+1;
[~,~,~,N] = RobotInit();
X_mont = cell(M,1);
X_hat_CD_mont = cell(M,1); X_hat_SD_mont = cell(M,1);
P_CD_mont = cell(M,1); P_SD_mont = cell(M,1);
W_mont=cell(M,1);
for Monte_index=1:M
    X_mont{Monte_index,1} = cell(k_f,N);
    X_hat_CD_mont{Monte_index,1} = cell(k_f,N);
    X_hat_SD_mont{Monte_index,1} = cell(k_f,N);
    P_CD_mont{Monte_index,1} = cell(k_f,1);
    P_SD_mont{Monte_index,1} = cell(k_f,1);
    for time_index = 1:k_f
        P_CD_mont{Monte_index,1}{time_index,1} = cell(N,N);
        P_SD_mont{Monte_index,1}{time_index,1} = cell(N,N);
    end
end
save .\AnalysisData.mat X_mont X_hat_CD_mont X_hat_SD_mont P_CD_mont P_SD_mont M;

%% Noise generation
[indep_noise{1:(N+N+5*N)*M}] = RandStream.create('mrg32k3a','NumStreams',7*N*M);
% [indep_noise{1:(N+N+5*N)*M}] = RandStream.create('mlfg6331_64','NumStreams',7*N*M);

% 7 types of random noise
Noise_propa_w_mont=cell(M,1);
Noise_propa_v_mont=cell(M,1);
Noise_extero_x_mont=cell(M,1);
Noise_extero_y_mont=cell(M,1);
Noise_extero_phi_mont=cell(M,1);
Noise_extero_ax_mont=cell(M,1);
Noise_extero_ay_mont=cell(M,1);
r1 = -0.5; r2 = 0.5;

for Monte_index=1:M
    for UID=1:N
        noise_index=UID+(Monte_index-1)*(N+N+5*N);
        Noise_propa_w_mont{Monte_index}(:,UID)=randn(indep_noise{noise_index},length_noise,1);
        Noise_propa_v_mont{Monte_index}(:,UID)=randn(indep_noise{noise_index+N},length_noise,1);
        Noise_extero_x_mont{Monte_index}(:,UID)=randn(indep_noise{noise_index+2*N},length_noise,1);
        Noise_extero_y_mont{Monte_index}(:,UID)=randn(indep_noise{noise_index+3*N},length_noise,1);
        Noise_extero_phi_mont{Monte_index}(:,UID)=randn(indep_noise{noise_index+4*N},length_noise,1);
        Noise_extero_ax_mont{Monte_index}(:,UID)=randn(indep_noise{noise_index+5*N},length_noise,1);
        Noise_extero_ay_mont{Monte_index}(:,UID)=randn(indep_noise{noise_index+6*N},length_noise,1);
    end
W_mont{Monte_index}=(r1 + (r2-r1).*rand(1,N));
end

save .\velocity_profile_mont.mat W_mont
save .\noise_profile_mont.mat  Noise_propa_v_mont   Noise_propa_w_mont ...
                        Noise_extero_x_mont  Noise_extero_y_mont Noise_extero_phi_mont ...
                        Noise_extero_ax_mont Noise_extero_ay_mont;
end