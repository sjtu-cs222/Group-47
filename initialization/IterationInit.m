%%%%%%%%%%Created by Qi Yan%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Last Revised August 2018%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [delta,k_f,w_m,w_c,lambda,RelMea_Table,UpdateOrder,opt_para] = IterationInit(missed_onoff)
% STATIC iteration parameters initialization, DYNAMIC parameters such as
% variance of noise are determined elsewhere.
% delta                 : stepsize
% k_f                   : total number of steps
% w_m                   : weights for mean calculation (UKF)
% w_c                   : weights for covarance calculation (UKF)
% RelMea_Table          : relative measurements occurrence table
% UpdateOrder           : update order for sequential updating

%% computation settings
delta = 0.1; % stepsize
t_f = 100; % total time
k_f = floor(t_f/delta); % total number of steps
% com_rate = 0.15;
% com_rate = 0.3;
opt_para = [2 2 2 25];
%% UKF weights settings
[~,~,~,N] = RobotInit();
n_x = 2*N; % number of all state variables

alpha = 1; kappa = (3-n_x); % determining parameters for UKF weights
beta = 2;
lambda = alpha^2*(n_x+kappa)-n_x;

w_m = zeros(1,2*n_x+1); w_c = zeros(1,2*n_x+1);
w_m(1) = lambda/(n_x+lambda);
w_c(1) = lambda/(n_x+lambda) + (1-alpha^2+beta);
for i=2:2*n_x+1
    w_m(i) = 1/(2*(n_x+lambda));
    w_c(i) = 1/(2*(n_x+lambda));
end

%% relative measurement times
% missed_onoff=0: Every one updates
% missed_onoff=1: Only two robots making relative measurements update (loosely coupled)
% missed_onoff=2: No CL

if nargin == 0
    missed_onoff = 0;
end

UpdateOrder = [1:N]; %this is the priority list for updating (sequential updating)

% each element(matrix) in RelMea_Table:
% [ mea. start iteration   , mea. end iteration   ;
%   identity of robot_a_1  , identity of robot_b_1;
%   identity of robot_a_2  , identity of robot_b_2;]

if missed_onoff == 0
%     %All the robots recieve the update message
%     RelMea_Table={[floor(10/delta)+1     floor(60/delta)  ; 2 1 ; 3 1; 4  1; 10 5; 11 6; 12 7; 13 7; 14 8; 15 9],...
%                   [floor(60/delta)+1     floor(110/delta) ; 5 2 ; 6 2; 7  3; 8  3; 9  4; 1 1],...
%                   [floor(110/delta)+1    floor(160/delta) ; 2 1 ; 3 1; 4  1; 10 2; 14 4; 15 5; 13 13],...
%                   [floor(160/delta)+1    floor(210/delta) ; 6 2 ; 7 3; 8  3; 12 1; 7  7],...
%                   [floor(210/delta)+1    floor(260/delta) ; 7 1 ;15 8; 1 13; 9  5],...
%                   [-1 -1;0 0]};
% RelMea_Table={
%     [-1 -1;0 0]};

% RelMea_Table={[floor(20/delta)+1 floor(22/delta);1 1],...
%     [-1 -1;0 0]};
RelMea_Table={[floor(20/delta)+1 floor(26/delta);1 2;2 3;3 1],...
    [floor(40/delta)+1 floor(46/delta);1 2;2 3;3 1],...
    [floor(60/delta)+1 floor(66/delta);1 2;2 3;3 1],...
    [floor(80/delta)+1 floor(86/delta);1 2;2 3;3 1],...
    [-1 -1;0 0]};
% RelMea_Table={[floor(50/delta)+1 floor(90/delta);3 3; 2 3; 3 2; 1 3],...
%     [-1 -1;0 0]}; % 7, 4, 5, 2
% RelMea_Table={[floor(50/delta)+1 floor(90/delta); 1 2; 1 3; 2 1; 2 3; 3 2; 3 1; 3 3],...
%     [-1 -1;0 0]}; % all in
elseif missed_onoff == 1
    %Only master and landmark robots update
    t_zero = zeros(1,N-1);
%     RelMea_Table={[1    1  t_zero; 3 4,setdiff(1:N,[3 4]) 0; 5 5,setdiff(1:N,[5 5])  ;2 1,setdiff(1:N,[2 1]) 0],...
%                   [10   12 t_zero; 3 4,setdiff(1:N,[3 4]) 0; 5 5,setdiff(1:N,[5 5])  ;2 1,setdiff(1:N,[2 1]) 0],...
%                   [50   53 t_zero; 2 3,setdiff(1:N,[2 3]) 0; 4 1,setdiff(1:N,[4 1]) 0;5 5,setdiff(1:N,[5 5])  ],...
%                   [61  100 t_zero; 4 3,setdiff(1:N,[4 3]) 0; 1 1,setdiff(1:N,[1 1])  ;2 5,setdiff(1:N,[2 5]) 0],...
%                   [110 145 t_zero; 5 3,setdiff(1:N,[5 3]) 0; 4 4,setdiff(1:N,[4 4])  ;1 2,setdiff(1:N,[1 2]) 0],...
%                   [150 250 t_zero; 4 4,setdiff(1:N,[4 4])  ; 5 3,setdiff(1:N,[5 3]) 0;2 2,setdiff(1:N,[2 2])  ],...
%                   [260 300 t_zero; 1 5,setdiff(1:N,[1 5]) 0; 3 3,setdiff(1:N,[3 3])  ;2 4,setdiff(1:N,[2 4]) 0],...
%                   [-1 -1;0 0]};
    RelMea_Table={[floor(10/delta)+1     floor(60/delta)  t_zero; 2 1 setdiff(1:N,[2 1]) 0 ; 5 1 setdiff(1:N,[5 1]) 0],...
                  [floor(60/delta)+1     floor(110/delta) t_zero; 3 2 setdiff(1:N,[3 2]) 0],...
                  [floor(110/delta)+1    floor(160/delta) t_zero; 1 1 setdiff(1:N,[1 1])],...
                  [floor(160/delta)+1    floor(210/delta) t_zero; 2 2 setdiff(1:N,[2 2])],...
                  [floor(210/delta)+1    floor(260/delta) t_zero; 3 1 setdiff(1:N,[3 1]) 0],...
                  [-1 -1;0 0]};
elseif missed_onoff == 2
    % No CL
    clear RelMea_Table;
    RelMea_Table={[-1 -1;0 0]};
end

end