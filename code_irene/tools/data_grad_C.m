function gradient_data = data_grad_C(A, C, Psi_fcts, ...
    lump_mass_mat_M, G, v, M, k, r, W, mu_3, mu_45)
eta = arrayfun(@(t) 1/2 * (tanh(2*t - 1) + 1), v);
B = Psi_fcts' * lump_mass_mat_M * (G .* eta); %[k, q]
H = C*A - B; %[k,q]
mask = eye(k); %[k,k]
yoda = @(p,q) sum(sum(H .* (mask(:,p) * A(q,:)), 1) ...
     ./ vecnorm(H), 2);
pq_comb = combvec(1:k, 1:k);
gradient_data = arrayfun(yoda, pq_comb(1,:), pq_comb(2,:)); %[k,k]
gradient_data = reshape(gradient_data, k, k); %[k,k]
gradient_data = norm()
end