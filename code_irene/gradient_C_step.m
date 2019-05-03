function gradient = gradient_C_step(M, G, A, C, v, W, r, par)

%% data term
eta = arrayfun(@(t) 1/2 * (tanh(2*t - 1) + 1), v);
B = M.evecs' * M.LMass * (G .* eta); %[k, q]
H = C*A - B; %[k,q]
mask = eye(par.k); %[k,k]
yoda = @(p,q) sum(sum(H .* (mask(:,p) * A(q,:)), 1) ...
     ./ vecnorm(H), 2);
pq_comb = combvec(1:par.k, 1:par.k);
gradient_data = arrayfun(yoda, pq_comb(1,:), pq_comb(2,:)); %[k,k]
gradient_data = reshape(gradient_data, par.k, par.k); %[k,k]


%% mu_3 term
d_mu_3_term = 2 * C .* (W.^2); %[k,k]
 
%% mu_45 term
D = cat(2, ones(1,r), zeros(1, par.k-r)); %[1,k]
d_mu_45_term = 4 * ((C*(C')*C) - D .* C); %[k,k]

%% gradient all together
gradient = gradient_data + par.mu_3 * d_mu_3_term ...
    + par.mu_45 * d_mu_45_term;
end 