%%
clc 
close all
clear all

% load a shape with a boundary
N = load_obj('device2-7.obj');
%N = load_off("sphere.off");
%load("wolf0_A0.40_H25.mat");
% N = rename_fields(N, "SGP");
k = 300;

%% così calcolo solo gli autovalori
% NOTA: in questo modo non posso ottenere l'operatore di schroedinger

[N.W,~,N.A] = calc_LB_FEM(N);
boundary = calc_boundary_edges(N.TRIV);
boundary = unique(boundary(:));
inside = setdiff(1:N.n, boundary);
[N.evecs1, N.evals1] = dirichlet_efcts(N, inside, k);
[N.evals1, idx] = sort(N.evals1);
N.evecs1 = N.evecs1(:,idx);


%% Laplaciano codice precedente

W2 = sparse(N.n, N.n);
l = size(boundary); l = l(1);
W2(boundary, boundary) = eye(length(boundary));
W2(inside,:) = N.W(inside, :);
N.W2 = W2;
% A2 = sparse(N.n, N.n);
% A2(inside,inside) = N.A(inside, inside);
% N.A2 = A2;
[N.evecs2, N.evals2] = eigs(N.W2, N.A, k, 'SM');
N.evals2 = diag(N.evals2);
[N.evals2, idx] = sort(N.evals2);

%% Laplaciano con bordo di Dirichlet debuggato

W3 = sparse(N.n, N.n);
l = size(boundary); l = l(1);
W3(boundary, boundary) = eye(length(boundary));
W3(inside, inside) = N.W(inside, inside);
N.W3 = W3;
A3 = sparse(N.n, N.n);
A3(inside,inside) = N.A(inside, inside);
N.A3 = A3;
[N.evecs3, N.evals3] = eigs(N.W3, N.A3, k, 1e-10);
[N.evals3, idx] = sort(diag(N.evals3));

%% li plotto entrambi
figure, trisurf(N.TRIV, N.VERT(:, 1), N.VERT(:, 2), N.VERT(:, 3));

figure, bar([N.evals1 N.evals2 N.evals3],...
    'LineWidth', 1)
legend(...
    'Dirichlet eigenvalues', ...
    'Laplacian with Dirichlet b.c. prima', ...
    'Laplacian with Dirichlet b.c. dopo', ...
    'Location', 'NorthWest')
xlabel('Index')
ylabel('Eigenvalue')