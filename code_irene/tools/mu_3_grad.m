function d_mu_3_term = mu_3_grad(A, C, Psi_fcts, ...
    lump_mass_mat_M, G, v, M, k, r, W, mu_3, mu_45)
    % mu_3 term
    d_mu_3_term = 2 * C .* (W.^2); %[k,k]

end