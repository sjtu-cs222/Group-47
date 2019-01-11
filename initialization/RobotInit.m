%%%%%%%%%%created by Qi Yan%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Last Revised August 2018%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X0,X_hat0,P0,N] = RobotInit()
%% Robots initialization
% X0     : initial true poses
% X_hat0 : initial estimated poses and covariances
% P0     : initial collective covariances
% N      : number of robots

%% Robot number
N = 10;
R2D = 180/pi;
D2R = pi/180;

%% True system pose X{ID}
X0 = cell(1,N);
X0{1} = [60,70]';
X0{2} = [40,40]';
X0{3} = [60,40]';
X0{4}=[80,40]';
X0{5}=[20,10]';
X0{6}=[40,10]';
X0{7}=[60,10]';
X0{8}=[80,10]';
X0{9}=[100,10]';
X0{10}=[20,-20]';
% X0{11}=[40,-20,-pi/2]';
% X0{12}=[60,-20,-pi/2]';
% X0{13}=[80,-20,-pi/2]';
% X0{14}=[100,-20,-pi/2]';
% X0{15}=[120,-20,-pi/2]';

%%% Estimated pose Xhat{ID}
X_hat0=cell(1,N);
for i = 1:N
    X_hat0{i}=X0{i}*(0e-3*randn+1);
end

%% Robot Covariance P{ID}
P_coll = cell(1,N);
for i = 1:N
%     P_coll{i} = diag([0.25,0.25,1*D2R]);
    P_coll{i} = 0.02^2*eye(2);
end
P0 = blkdiag(P_coll{:});

end