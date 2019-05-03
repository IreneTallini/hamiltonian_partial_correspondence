function c = cost_v_step(N, M, G, A, C, v, par)

%% Compute B = <M.evecs, G .* eta(v)>_M:
% it's necessary to recompute it at each step as it depends on v
eta = arrayfun(@(t) 1/2 * (tanh(2*t - 1) + 1), v); % soft membership function 
                                                % (rescales v from (-inf, inf) to (0, 1))
B = M.evecs' * M.LMass * (G .* eta); %[k, q]

%% Data term: L1,2 norm of CA - B
data_term = sum(vecnorm(C*A - B));

%% Regularization constants (see def (13) of paper [1])
%% mu_1 term: promotes similarity between area of N and area of M masked via eta(v).
area_N = sum(diag(N.LMass));
int_eta = sum(M.LMass * eta);
mu_1_term = (area_N - int_eta)^2;

%% mu_2 term: promotes short boundary of M masked via eta(v).
E = vecnorm(M.VERT(M.TRIV(:,2),:) - M.VERT(M.TRIV(:,1),:), 2, 2).^2; 
F = dot(M.VERT(M.TRIV(:,2),:) - M.VERT(M.TRIV(:,1),:), M.VERT(M.TRIV(:,3),:) - M.VERT(M.TRIV(:,1),:), 2);
G = vecnorm(M.VERT(M.TRIV(:,3),:) - M.VERT(M.TRIV(:,1),:), 2, 2).^2;
v_alpha = v(M.TRIV(:,2)) - v(M.TRIV(:,1));
v_beta = v(M.TRIV(:,3)) - v(M.TRIV(:,1));
D = (v_alpha.^2 .* G - 2 * v_alpha .* v_beta .* F + ...
    v_beta.^2 .* E).^(1/2) ; 
xi = @(t) exp(-tanh(2*t - 1)/(4 * par.sig^2));
xi_sum = arrayfun(xi, v(M.TRIV(:,1))) + ...
    arrayfun(xi, v(M.TRIV(:,2))) + ...
    arrayfun(xi, v(M.TRIV(:,3)));
mu_2_term = 1/6 * sum(D .* xi_sum);

%% cost all together
c = data_term + par.mu_1 * mu_1_term + par.mu_2 * mu_2_term;
end