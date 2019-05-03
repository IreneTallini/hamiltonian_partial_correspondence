function gradient = gradient_v_step(N, M, G, A, C, v, par)

eta = arrayfun(@(t) 1/2 * (tanh(2*t - 1) + 1), v);
B = M.evecs' * M.LMass * (G .* eta); %[k, q]

%% data term 
H = C*A - B; %[k,q]
d_eta = arrayfun(@(t) (sech(2*t - 1))^2, v);
Psi_fcts_T = M.evecs'; %[k, M.n]
d = full(diag(M.LMass));
yoda = @(p) sum(sum(H .* (Psi_fcts_T(:,p) *...
    ((-1) * d(p) * d_eta(p) * G(p,:))), 1) ./ vecnorm(H), 2); 
gradient_data = arrayfun(yoda, 1:M.n)';

%% mu_1 term
area_N = sum(diag(N.LMass)); %[1]
eta = arrayfun(@(t) 1/2 * (tanh(2*t - 1) + 1), v);
int_eta = sum(M.LMass * eta); %[1]
d_int_eta = diag(M.LMass) .* d_eta; 
gradient_mu1_term = (-2) * (area_N - int_eta) * d_int_eta; %[n,1]

%% mu_2 term
xi = @(t) exp(-tanh(2*t - 1)/(4 * par.sig^2)); 
d_xi = arrayfun(@(t) ...
    -(((sech(2*t - 1)^2)/(2 * par.sig^2)) * ...
    exp(-tanh(2*t - 1)/(4 * par.sig^2))), v);
E = vecnorm(M.VERT(M.TRIV(:,2),:) - M.VERT(M.TRIV(:,1),:), 2, 2).^2; 
F = dot(M.VERT(M.TRIV(:,2),:) - M.VERT(M.TRIV(:,1),:), ...
    M.VERT(M.TRIV(:,3),:) - M.VERT(M.TRIV(:,1),:), 2);
G = vecnorm(M.VERT(M.TRIV(:,3),:) - M.VERT(M.TRIV(:,1),:), 2, 2).^2;
v_alpha = v(M.TRIV(:,2)) - v(M.TRIV(:,1));
v_beta = v(M.TRIV(:,3)) - v(M.TRIV(:,1));
D = (v_alpha.^2 .* G - 2 * v_alpha .* v_beta .* F ...
    + v_beta.^2 .* E) .^(1/2); %[m]
vert_to_tri_1 = sparse(1:M.m, M.TRIV(:,1), ones(1,M.m), M.m, M.n)';
vert_to_tri_2 = sparse(1:M.m, M.TRIV(:,2), ones(1,M.m), M.m, M.n)';
vert_to_tri_3 = sparse(1:M.m, M.TRIV(:,3), ones(1,M.m), M.m, M.n)';
K_1 = D.^(-1) .* ...
    ((v(M.TRIV(:,1)) - v(M.TRIV(:,2))) .* (G - F) + ...
    (v(M.TRIV(:,1)) - v(M.TRIV(:,3))) .* (E - F));
K_2 = D.^(-1) .* ...
    ((v(M.TRIV(:,2)) - v(M.TRIV(:,1))) .* G - ...
    (v(M.TRIV(:,3)) - v(M.TRIV(:,1))) .* F);
K_3 = D.^(-1) .* ...
    ((v(M.TRIV(:,3)) - v(M.TRIV(:,1))) .* E - ...
    (v(M.TRIV(:,2)) - v(M.TRIV(:,1))) .* F);  
K = [K_1 K_2 K_3]; %[m,3] la prima colonna contiene K rispetto al primo vertice etc...
xi_sum = arrayfun(xi, v(M.TRIV(:,1))) + ...
    arrayfun(xi, v(M.TRIV(:,2))) + ...
    arrayfun(xi, v(M.TRIV(:,3))); %[1,m]
gradient_mu2_term = 1/6 * ...
    (vert_to_tri_1 * (xi_sum .* K(:,1) + D .* d_xi(M.TRIV(:,1))) + ...
    vert_to_tri_2 * (xi_sum .* K(:,2) + D .* d_xi(M.TRIV(:,2))) + ...
    vert_to_tri_3 * (xi_sum .* K(:,3) + D .* d_xi(M.TRIV(:,3))));
gradient_mu2_term(isinf(gradient_mu2_term) | isnan(gradient_mu2_term)) = 0;

%% gradient all together
gradient = gradient_data + par.mu_1 * gradient_mu1_term + par.mu_2 * gradient_mu2_term;

end 