clear all; close all; clc;

%% load full meshes from TOSCA
mesh_full_path = "centaur0";
mesh_partial_path = "centaur1";
mesh_full = load(strcat(mesh_full_path, '.mat'));
mesh_partial = load(strcat(mesh_partial_path, '.mat'));

%% change from TOSCA nomenclature to mine
M = rename_fields(mesh_full.surface, "TOSCA");
N_tmp = rename_fields(mesh_partial.surface, "TOSCA");

%% cut N with plane. this will be the partial shape
[N, map_N_to_M] = cut_with_plane([0 0 0], [0 -30 10], N_tmp);

%% parameters
par.k = 100; % number of eigenvalues/eigenfunctions
par.sig = 0.03; % see "Weight matrix (?3-term)" section of paper [1]
par.mu_1 = 1; 
par.mu_2 = 1e+2;
par.mu_3 = 1;
par.mu_45 = 1e+3;
par.num_it = 5; % number of outer iterations of the alternating scheme

%% compute SHOT descriptors 
F = calc_shot(N.VERT', N.TRIV', 1:N.n, 9, 10, 3)'; %[N.n, 320]
F(:, ~any(F,1)) = []; % Delete all zero descriptors
G = calc_shot(M.VERT', M.TRIV', 1:M.n, 9, 10, 3)'; %[M.n, 320]
G(:, ~any(G,1)) = [];  % Delete all zero descriptors

%% shape matching
[M, N, v, C] = partial_shape_matching(M, N, F, G, par);

%% print results
print_results;