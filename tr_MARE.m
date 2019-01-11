function tr = tr_MARE(Pi_init,Q,H,W)
pi_cell = cell(10,1); Tr_pi = zeros(10,1);
pi_cell{1} = Pi_init;
dim = length(H);
for i=1:1e6
    pi_cell{i+1} = pi_cell{i} - pi_cell{i}*H'*(W.*(H*pi_cell{i}*H' + eye(dim)))^-1*H*pi_cell{i};
    pi_cell{i+1} = pi_cell{i+1} + Q;
    Tr_pi(i) = trace(pi_cell{i});
    if i > 100
        if abs((Tr_pi(i) - Tr_pi(i-1))/Tr_pi(i)) < 1e-3
            break;
        end
    end
end
tr = Tr_pi(end);
end