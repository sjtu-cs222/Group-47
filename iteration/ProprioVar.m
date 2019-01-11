%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%created by Solmaz S. Kia%%%%%%%%%%%%%%%
%%%%%%%%%%Last Revised April 2014%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%UCSD%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sigma_V,sigma_phi] = ProprioVar(V,k)
%% uncertatinty for true robot motion (proprioceptive sensor)
persistent dt N W;
if isempty(dt) && isempty(N)
    [dt,~] = IterationInit();
%     [~,~,~,N] = RobotInit();
    N = length(V);
    [~,W] = TrajGen(0);
end
%%
% linear velocity noise (standard deviation)
% sigma_V(1)=abs(0.1*V(1));%sqrt(0.04);
% sigma_V(2)=abs(0.05*V(2));%sqrt(0.02);
% sigma_V(3)=abs(0.1*V(3));%sqrt(0.01);
% sigma_V(4)=dt*abs(0.05*V(4));%sqrt(0.05);
% sigma_V(5)=dt*abs(0.05*V(5));%sqrt(0.04);
sigma_V = zeros(1,N);
for i = 1:N
    sigma_V(i)=abs(sqrt(5.075)*V(i));
end

% angular velocity noise (standard deviation)
% sigma_phi(1)=(0.5*pi/180);%sqrt(0.001);
% sigma_phi(2)=(0.5*pi/180);%sqrt(0.01);
% sigma_phi(3)=(0.5*pi/180);%sqrt(0.01);
% sigma_phi(4)=(0.5*pi/180);%sqrt(0.01);
% sigma_phi(5)=0.5*pi/180;%sqrt(0.02);
sigma_phi = zeros(1,N);
for i = 1:N
    %     sigma_phi(i) = (3*pi/180); % 3 degree
    %     sigma_phi(i) = 10*pi/180 + abs(W(i)*0.35); % 3 degree
    sigma_phi(i) = sqrt(0.345);
end
end