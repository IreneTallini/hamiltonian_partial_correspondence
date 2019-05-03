function gradient_mu1_term = mu1_grad(...
    A, C, Psi_fcts, G, M, v, lump_mass_mat_M, lump_mass_mat_N, mu_1, mu_2, sig)
%mu_1 term
area_N = sum(diag(lump_mass_mat_N)); %[1]
eta = arrayfun(@(t) 1/2 * (tanh(2*t - 1) + 1), v);
int_eta = sum(lump_mass_mat_M * eta); %[1]
d_int_eta = diag(lump_mass_mat_M) .* d_eta; 
gradient_mu1_term = (-2) * (area_N - int_eta) * d_int_eta; %[n,1]
end