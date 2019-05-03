function mu_1_term = mu_1_term_v(C, A, Psi_fcts, ...
    lump_mass_mat_M, lump_mass_mat_N, G, v, M, mu_1, mu_2, sig)
area_N = sum(diag(lump_mass_mat_N));
eta = arrayfun(@(t) 1/2 * (tanh(2*t - 1) + 1), v);
int_eta = sum(lump_mass_mat_M * eta);
mu_1_term = (area_N - int_eta)^2;
end