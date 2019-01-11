clear;
close all
F = 0.1*blkdiag(2,3);
Q = @(t) blkdiag(1,sin(t*pi/2));
P0 = 1*eye(2); DP = 2*eye(2);
P_init = P0 + DP;
P = cell(2,3); Phi = cell(1,2);
P{1,1} = P_init; P{1,2} = P0; P{1,3} = DP; Phi{1} = eye(2);
for i = 1:10
    P{i+1,1} = F*P{i,1}*F' + Q(i);
%     P{i+1,1} = F*P{i,1}*F';
    P{i+1,2} = F*P{i,2}*F' + Q(i);
%     P{i+1,2} = F*P{i,2}*F';
    P{i+1,3} = F*P{i,3}*F';
    Phi{i+1} = F*Phi{i};
end
% error = P{end,1} - P{end,2} - Phi{end}*DP*Phi{end}'
error = P{end,1} - P{end,2} - P{end,3}