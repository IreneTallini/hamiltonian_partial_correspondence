%%
clc 
close all
clear all

% load a shape with a boundary
N = load_off('./cuts_david_shape_1.off');
k = 100;

%% così calcolo solo gli autovalori
% NOTA: in questo modo non posso ottenere l'operatore di schroedinger

[N.W,~,N.A] = calc_LB_FEM(N);
boundary = calc_boundary_edges(N.TRIV);
boundary = unique(boundary(:));
inside = setdiff(1:N.n, boundary);
[N.evecs1, N.evals1] = dirichlet_efcts(N, inside, k);
[N.evals1, idx] = sort(N.evals1);
N.evecs1 = N.evecs1(:,idx);

%% Laplaciano con bordo di Dirichlet debuggato

W2 = sparse(N.n, N.n);
l = size(boundary); l = l(1);
W2(boundary, boundary) = eye(length(boundary));
% Qui le colonne corrispondenti a vertici di bordo vanno lasciate a zero
W2(inside,inside) = N.W(inside, inside);
N.W2 = W2;
A2 = sparse(N.n, N.n);
A2(inside,inside) = N.A(inside, inside);
N.A2 = A2;
[N.evecs2, N.evals2] = eigs(N.W2, N.A2, k, 1e-10);
[N.evals2, idx] = sort(diag(N.evals2));

%% li plotto entrambi

figure, bar([...
    N.evals1 ...
    N.evals2],...
    'LineWidth', 1)
legend(...
    'Dirichlet eigenvalues', ...
    'Laplacian with Dirichlet b.c.', ...
    'Location', 'NorthWest')
xlabel('Index')
ylabel('Eigenvalue')