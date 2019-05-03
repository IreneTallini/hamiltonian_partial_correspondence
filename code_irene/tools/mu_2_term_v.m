function mu_2_term = mu_2_term_v(C, A, Psi_fcts, ...
    lump_mass_mat_M, lump_mass_mat_N, G, v, M, mu_1, mu_2, sig)
E = vecnorm(M.VERT(M.TRIV(:,2),:) - M.VERT(M.TRIV(:,1),:), 2, 2).^2; 
F = dot(M.VERT(M.TRIV(:,2),:) - M.VERT(M.TRIV(:,1),:), M.VERT(M.TRIV(:,3),:) - M.VERT(M.TRIV(:,1),:), 2);
G = vecnorm(M.VERT(M.TRIV(:,3),:) - M.VERT(M.TRIV(:,1),:), 2, 2).^2;
v_alpha = v(M.TRIV(:,2)) - v(M.TRIV(:,1));
v_beta = v(M.TRIV(:,3)) - v(M.TRIV(:,1));
D = (v_alpha.^2 .* G - 2 * v_alpha .* v_beta .* F + ...
    v_beta.^2 .* E).^(1/2) ; 
xi = @(t) exp(-tanh(2*t - 1)/(4 * sig^2));
xi_sum = arrayfun(xi, v(M.TRIV(:,1))) + ...
    arrayfun(xi, v(M.TRIV(:,2))) + ...
    arrayfun(xi, v(M.TRIV(:,3)));
mu_2_term = 1/6 * sum(D .* xi_sum);
end