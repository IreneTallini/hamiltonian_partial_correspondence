function data_term = data_term_v(C, A, Psi_fcts, ...
    lump_mass_mat_M, lump_mass_mat_N, G, v, M, mu_1, mu_2, sig)
eta = arrayfun(@(t) 1/2 * (tanh(2*t - 1) + 1), v);
B = Psi_fcts' * lump_mass_mat_M * (G .* eta); %[k, q]
%Data term
data_term = sum(vecnorm(C*A - B));
end