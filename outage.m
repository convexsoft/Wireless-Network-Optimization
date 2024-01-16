function [p,power_evolution]=outage(G,n,beta,pmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem parameter initialization
% L: number of users.
% G: channel gain.
% F: nonnegative matrix. F_{lj} = G_{lj} if l ~= j, and F_{lj} = 0 if l = j
% n: noise vector
% v: nonnegative vector. v_l = n_l/G_{ll}
% beta: a vector that is a weight assigned to links to reflect priority.
% pmax: upper bound of the total power constraints.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    F = diag(diag(1./G))*(G-diag(diag(G)));
    v = n./diag(G); 
    [L,~] = size(G);
    p = rand(L,1);
    power_evolution = [p'];
    tol = 10e-7;
    err = 1;
   
    while err > tol
        p_temp = p;
        p = (v.*beta./p + log(prod(1+diag(beta./p)*F*diag(p),2))).*p;
        power_evolution = [power_evolution;(p/sum(p)*pmax)'];
        err = norm(p-p_temp,2);
    end
    p = p/sum(p)*pmax;
end

