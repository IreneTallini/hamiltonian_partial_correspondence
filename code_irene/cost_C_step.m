function c = cost_C_step(M, G, A, C, v, W, r, par)

%% Compute B = <M.evecs, G .* eta(v)>_M:
% it's necessary to recompute it at each step as it depends on v
eta = arrayfun(@(t) 1/2 * (tanh(2*t - 1) + 1), v);
B = M.evecs' * M.LMass* (G .* eta); %[k, q]

%% Data term: L1,2 norm of CA - B
data_term = sum(vecnorm(C*A - B));

%% Regularization constants (see def (14) of paper [1])
%% mu_3 term: promotes slanted diagonal shape of C
mu_3_term = norm(C .* W, 'fro')^2;

%% mu_4,5 term: promotes orthogonality of a submatrix of C'C. 
% NON HO IDEA DEL SENSO DI QUESTA COSA.
off_diagonal_penalty = norm(C'*C, 'fro')^2 - trace((C'*C).^2);
D = cat(2, ones(1,r), zeros(1, par.k-r)); %[1,k]
diagonal_penalty = trace((C'*C - D).^2);
mu_45_term = off_diagonal_penalty + diagonal_penalty;

%% cost all together
c = data_term + par.mu_3 * mu_3_term + par.mu_45 * mu_45_term ;

end