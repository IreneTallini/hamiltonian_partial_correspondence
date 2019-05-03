function mu_45_term = mu_45_term_C(C, A, Psi_fcts, lump_mass_mat_M,...
    G, v, r, k, mu_3, mu_45, W)
off_diagonal_penalty = norm(C'*C, 'fro')^2 - trace((C'*C).^2);
D = cat(2, ones(1,r), zeros(1, k-r)); %[1,k]
diagonal_penalty = trace((C'*C - D).^2);
mu_45_term = off_diagonal_penalty + diagonal_penalty;
end