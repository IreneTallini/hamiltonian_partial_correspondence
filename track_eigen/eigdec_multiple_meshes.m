function meshes = eigdec_multiple_meshes(meshes, n_eigen)
% EIGED_MULTIPLE_MESHES Computes eigenvalues and eigenfunction of a list of
% meshes

s = size(meshes); s = s(2);
    for i = 1:s
        [meshes(i).Stiff, meshes(i).Mass, meshes(i).LMass] = calc_LB_FEM(meshes(i));
        [meshes(i).evecs, meshes(i).evals] = eigs(meshes(i).Stiff, meshes(i).LMass, n_eigen, -1e-5);
    end
end