function d_mu_45_term = mu_45_grad (A, C, Psi_fcts, ...
    lump_mass_mat_M, G, v, M, k, r, W, mu_3, mu_45)
D = cat(2, ones(1,r), zeros(1, k-r)); %[1,k]
d_mu_45_term = 4 * ((C*(C')*C) - D .* C); %[k,k]
end