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
inside = calc_boundary_edges(N.TRIV);
inside = setdiff(1:N.n, unique(inside(:)));
[N.evecs1, N.evals1] = dirichlet_efcts(N, inside, k);
[N.evals1, idx] = sort(N.evals1);
N.evecs1 = N.evecs1(:,idx);

%% così calcolo il Laplaciano con bordo di dirichlet
% NOTA: questo procedimento potrebbe avere dei bug, a seconda di N.n e di k
% possono venire autovalori negativi. Se non sbaglia vengono uguali al caso
% precedente

N_obj = TriangleMesh(N);
N_obj.boundary_mode = 'dirichlet';
N.A_dir = N_obj.mass('barycentric');
N.W_dir = -N_obj.stiffness('dirichlet_full');

[N.evecs2, N.evals2] = eigs(N.W_dir, N.A, k, 'SM');
N.evals2 = diag(N.evals2);
[N.evals2, idx] = sort(N.evals2);
N.evecs2 = N.evecs2(:,idx);

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