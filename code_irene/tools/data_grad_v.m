function gradient_data = data_grad_v(...
    A, C, Psi_fcts, G, M, v, lump_mass_mat_M, lump_mass_mat_N, mu_1, mu_2, sig)
eta = arrayfun(@(t) 1/2 * (tanh(2*t - 1) + 1), v);
B = Psi_fcts' * lump_mass_mat_M * (G .* eta); %[k, q]

% data term 
H = C*A - B; %[k,q]
% d_eta = arrayfun(@(t) 1 - (tanh(2*t - 1))^2, v); %[M.n,1]
d_eta = arrayfun(@(t) (sech(2*t - 1))^2, v);
Psi_fcts_T = Psi_fcts'; %[k, M.n]
d = full(diag(lump_mass_mat_M));
yoda = @(p) sum(sum(H .* (Psi_fcts_T(:,p) *...
    ((-1) * d(p) * d_eta(p) * G(p,:))), 1) ./ vecnorm(H), 2); 
gradient_data = arrayfun(yoda, 1:M.n)';
end