q = 320;
S = fps_euclidean(N.VERT, q, 1);
F_del = eye(N.n); 
F_del = F_del(:, S); %[N.n, q]
G_del = eye(M.n);
G_del = G_del(:, map_N_to_M(S)); %[N.n, q]
eta = arrayfun(@(t) 1/2 * (tanh(2*t - 1) + 1), v);
A = N.evecs' * N.LMass* F_del; %[k, q]
B = M.evecs' * M.LMass * (G_del .* eta); %[k, q]
F_approx = N.evecs* A; %[M.n, q]
G_from_F_approx = M.evecs* C * A; %[M.n, q]


figure; 
trisurf(N.TRIV, N.VERT(:, 1), N.VERT(:, 2), N.VERT(:, 3), ...
    F_approx(:, 320));
shading interp;
title('F approx');
figure; 
trisurf(M.TRIV, M.VERT(:, 1), M.VERT(:, 2), M.VERT(:, 3), ...
    G_from_F_approx(:, 320));
shading interp;
title('G from F approx');
figure; 
trisurf(M.TRIV, M.VERT(:, 1), M.VERT(:, 2), M.VERT(:, 3), ...
    eta);
colormap hot;
shading interp;
title('mask');
colorbar;
figure('Name', "C matrix"); 
imagesc(C); colorbar;