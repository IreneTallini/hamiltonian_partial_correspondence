function mu_3_term = mu_3_term_C(C, A, Psi_fcts, lump_mass_mat_M,...
    G, v, r, k, mu_3, mu_45, W)
mu_3_term = norm(C .* W, 'fro')^2;
end