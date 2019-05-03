function [M, N, v, C] = partial_shape_matching(M, N, F, G, par)
% Implementation of paper "partial Functional Correspondance from Rodolà et
% al." [1].
%% compute Laplace eigenfunctions 
[N.Stiff, N.Mass, N.LMass] = calc_LB_FEM(N); %Computes cotangent FEM COME NELLE SLIDE DELLE LEZIONI
                                                % (QUINDI CON LE CONDIZIONI DI NEUMANN??)
                                                % IL CODICE è QUELLO CHE HA
                                                % DATO RODOLà A LEZIONE.
                                            
[M.Stiff, M.Mass, M.LMass] = calc_LB_FEM(M);
i = find(diag(full(M.LMass)) == 0); 
M.LMass(sub2ind(size(M.LMass), i, i)) = 1e-5; % Substitute zero elements in M.LMass with 1e-5.
                                               % NON SO BENE PERCHE SUCCEDA
                                               % STA COSA.
[N.evecs, N.evals] = eigs(N.Stiff, N.LMass, par.k, -1e-5); %[[N.n,k] [k, k]]
[M.evecs, M.evals] = eigs(M.Stiff, M.LMass, par.k, -1e-5); %[[M.n,k] [k, k]]

%% compute A coefficients of F
A = N.evecs' * N.LMass * F; %[k, q]

%% compute r: slanted diagonal inclination
% Essenzialmente l'unico momento in cui si usa r, ovvero il fatto che la
% matrice C è slanted diagonal, è per inizializzare C e e per calcolare la
% costante di regolarizzazione mu_3_term (vedi cost_C_step). Però se levi
% quell'informazione non funziona più.
max_M = max(diag(M.evals));
idx = find((diag(N.evals) < max_M) == 1);
r = idx(end);

%% compute W: matrix with zeros on the slanted diagonal and large values outside
W = W_matrix(par, r); %[k,k]

%% alternating scheme: first we optimize over C, then over v
% initialize C to W
C = W;
% initialize v to all ones
v = ones(M.n, 1); %[M.n,1] 

for i = 1:par.num_it
    % define problem for C optimization
    option.maxiter = 10000; 
    problem_C.M = euclideanfactory(par.k, par.k);
    problem_C.cost = @(C) cost_C_step(M, G, A, C, v, W, r, par);
    problem_C.egrad = @(C) gradient_C_step(M, G, A, C, v, W, r, par);
    % check if C gradient correct
    if i == 1
        checkgradient(problem_C);
    end
    % optimize C
    [C, ~, ~, ~] = conjugategradient(problem_C, C, option);
    
    % define problem for v optimization
    % NB: v can assume real values in (-inf, inf)
    problem_v.M = euclideanfactory(M.n, 1);
    problem_v.cost = @(v) cost_v_step(N, M, G, A, C, v, par);
    problem_v.egrad = @(v) gradient_v_step(N, M, G, A, C, v, par);
    
    % check if v gradient correct
    if i == 1
        checkgradient(problem_v);
    end
    % optimize v
    [v, ~, ~, ~] = conjugategradient(problem_v, v, option);
    
end

end