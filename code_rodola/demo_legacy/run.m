%
% This code accompanies the paper:
%
% "Partial Functional Correspondence"
% E.Rodola, L.Cosmo, M.M.Bronstein, A.Torsello, D.Cremers
% Computer Graphics Forum
%
% Please cite the paper above if you use this code.
%
% Written by L.Cosmo and E.Rodola
% Universita' Ca' Foscari Venezia, Italy
% TU Munich, Germany
% (c) 2015
%

clc, close all, clear all
addpath('./tools/')
addpath('./partial/')
addpath(genpath('./manopt/manopt/'))

% *************************************************************************
% Note that this demo code is temporary and is not optimized for speed. As
% such, it may take several minutes to converge. You can reduce 
% options.maxiter inside optimize_C.m and especially in optimize_v.m to 
% achieve some speed-up.
%
% A final version of the code will be made available soon.
% *************************************************************************

options = {};

options.n_eigen = 30;           % no. eigenpairs used
options.tv_sigma = 0.2;         % variance of indicator func. in TV term
options.tv_mean = 0.5;          % mean =
options.descr_max_radius = 5;   % max support radius for rigid descriptors (in % of diameter)
options.intinv_volstep = 0.3;   % voxelization step for integral invariants (in % of diameter)
options.shot_bins = 5;          % no. bins for SHOT descriptor
options.shot_min_neighs = 3;    % min no. neighbors for SHOT descriptor

% -------------------------------------------------------------------------
% Load the shapes
% -------------------------------------------------------------------------

M = load_off('./cat5.off');
M.S_tri = calc_tri_areas(M);

load('./cat0.mat');
N.S_tri = calc_tri_areas(N);

[M.evecs, M.evals, M.S, M.W] = calc_LB(M, options.n_eigen);
[N.evecs, N.evals, N.S, N.W] = calc_LB(N, options.n_eigen);

% -------------------------------------------------------------------------
% Compute dense descriptors
% -------------------------------------------------------------------------

diam = sqrt(sum(M.S_tri));
[F1, F2] = calc_dense_descriptors(N, diam, options);
[G1, G2] = calc_dense_descriptors(M, diam, options);
F = [F1 F2];
G = [G1 G2];
clear F1 F2 G1 G2

% -------------------------------------------------------------------------
% Match the partial shape N to the full shape M. This code reproduces the
% results shown in Figure 9 of the paper.
% -------------------------------------------------------------------------

v_init = ones(M.n,1);
C_init = [];

mu1 = 1;   % slanted-diagonal mask
mu2 = 1e3; % sub-orthogonality

figure
colors = create_colormap(M,M);
colormap(colors), plot_scalar_map(M, 1:M.n), axis off

for i=1:5
    
    [C, est_rank] = optimize_C(M, N, G, F, v_init, C_init, mu1, mu2);
    
    fprintf('Estimated rank: %d\n', est_rank);
    
%     figure
%     bar(sqrt(sum( (C*N.evecs'*N.S*F-M.evecs'*M.S*G).^2 )))
%     set(gca, 'xlim', [1 200])
%     title('Column-wise 2-norm (data term)')
    
    figure
    subplot(1,2,1), imagesc(C), axis equal; colorbar, title('C'); axis tight
    subplot(1,2,2), imagesc(C'*C), axis equal; colorbar, title('C^T C'); axis tight
    
    [~, matches] = run_icp(M, N, est_rank, C, 0);
    
    figure
    N2 = M; N2.VERT = M.VERT(matches,:);
    colors = create_colormap(N2,M);
    colormap(colors), plot_scalar_map(N, 1:N.n), axis off
    title(sprintf('Iteration %d\nCorrespondence before ICP',i))
    
    [Co, matches] = run_icp(M, N, est_rank, C, 1);
    
    figure
    N2 = M; N2.VERT = M.VERT(matches,:);
    colors = create_colormap(N2,M);
    colormap(colors), plot_scalar_map(N, 1:N.n), axis off
    title(sprintf('Iteration %d\nCorrespondence after ICP',i))
    
    C_init = [Co zeros(size(C,1),size(C,1)-size(Co,2))];
    
    v = optimize_v(M, N, G, F, C_init, 1, 1e2, options);
    v = (0.5*tanh(6*(v-0.5))+0.5);
    
    figure
    plot_scalar_map(M,v), title(sprintf('Iteration %d\nMask area agreement: %.2e',i,full(abs(sum(diag(N.S)) - sum(v.*diag(M.S))))))
    
    v_init = v;
    
end
