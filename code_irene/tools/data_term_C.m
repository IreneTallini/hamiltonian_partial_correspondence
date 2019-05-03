function data_term = data_term_C(C, A, Psi_fcts, lump_mass_mat_M,...
    G, v, r, k, mu_3, mu_45, W)
eta = arrayfun(@(t) 1/2 * (tanh(2*t - 1) + 1), v);
B = Psi_fcts' * lump_mass_mat_M * (G .* eta); %[k, q]
%Data term (assume B = B(nabla(v)))
data_term = sum(vecnorm(C*A - B));

end