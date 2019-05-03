function plot_functional_correspondance(M, N, C, indices, F)
% Input:
%   - M, N: meshes with fields for eigenfunctions, eigenvalues, stiffness and
%       mass matrix
%   - C: functional correspondence matrix
%   - indices: indices of the delta functions to plot 

% Compute coefficients of delta functions
del_M_coeff = M.evecs' * M.LMass * F; %[k q]
del_N_coeff = C * del_M_coeff; %[k q]
del_M = M.evecs * del_M_coeff;
del_N = N.evecs * del_N_coeff;

for i = indices
    figure;
    subplot(1,2,1); trisurf(M.TRIV, M.VERT(:,1), M.VERT(:,2), M.VERT(:,3), del_M(:, i)), ...
        axis equal; colorbar, title(strcat('delta of vertex ', num2str(i))); 
        axis tight; 
    subplot(1,2,2), trisurf(N.TRIV, N.VERT(:,1), N.VERT(:,2), N.VERT(:,3), del_N(:, i)),...
        axis equal; colorbar;
        title(strcat('delta of vertex ', num2str(i))); 
        axis tight;
end