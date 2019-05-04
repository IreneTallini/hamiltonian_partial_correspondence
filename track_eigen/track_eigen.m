function track_eigen(Ms, evec_ids)
% TRACK_EIGEN Plots tracking of a eigenvector along meshes.
%   TRACK_EIGEN(mesh_sequence, evecs_ids) tracks eigenvectors in evec_ids along 
%   the meshes in mesh_sequence.
%   
%   See also PLOT_FUNCTIONAL_CORRESPONDANCE

% Change evec_ids into a column vector
st = size(evec_ids);
if st(1) == 1
    evec_ids = evec_ids';
end
sm = size(Ms); sm = sm(2);
st = size(evec_ids); st = st(1);
track = zeros(st, sm, 'int8'); track(:,1) = evec_ids;

% Track eigenfunctions and plot C matrices
figure;
for i = 1 : sm - 1
    M = Ms(i);
    N = Ms(i+1);
    gt = [(1:1:M.n)' (1:1:M.n)'];
    TPhi = M.evecs(gt(:,2), :);
    C = TPhi' * N.LMass * N.evecs; %[n_eigen n_eigen]
    subplot(ceil(sqrt(sm)), ceil(sqrt(sm)), i);
    imagesc(C); colorbar; title(strcat("M = ", num2str(i), " vs N = ", num2str(i+1)));
    [~, tmp] = max(abs(C), [], 2);
    track(:,i+1) = tmp(track(:,i));
end

% Line plots of eigenfunctions along meshes
f1 = figure('Name', 'Tracking of eigenfunctions');
guidata(f1, 1);
plot(1:sm, track(1,:));
xlabel('Mesh'); ylabel(strcat("Eigenfunction ", num2str(1)));
annotation("arrow", [0.05 0.02], [0.5 0.5], "ButtonDownFcn", {@Scroll_plot, f1, sm, track, 'b'});
annotation("arrow", [0.95 0.98], [0.5 0.5], "ButtonDownFcn", {@Scroll_plot, f1, sm, track, 'f'});

% trisurf plot of eigenfunctions along meshes
M = Ms(1);
f2 = figure('Name', strcat("Eigenvector ", num2str(evec_ids(1)), " Mesh ", num2str(evec_ids(1))));
guidata(f2, [1 1]);
trisurf(M.TRIV, M.VERT(:,1), M.VERT(:,2), M.VERT(:,3), M.evecs(:,track(1,1))); 
axis equal;
annotation("arrow", [0.05 0.02], [0.5 0.5], "ButtonDownFcn", {@Scroll_mesh, f2, sm, st, Ms, track, 'l'});
annotation("arrow", [0.95 0.98], [0.5 0.5], "ButtonDownFcn", {@Scroll_mesh, f2, sm, st, Ms, track, 'r'});
annotation("arrow", [0.5 0.5], [0.95 0.98], "ButtonDownFcn", {@Scroll_mesh, f2, sm, st, Ms, track, 'u'});
annotation("arrow", [0.5 0.5], [0.05 0.02], "ButtonDownFcn", {@Scroll_mesh, f2, sm, st, Ms, track, 'd'});
end

% Callback function for trisurf plot
function Scroll_mesh(~, ~, f, sm, st, Ms, track, bf)
    % i,ud = evec
    % j,rl = mesh
    gd = guidata(f);
    i = gd(1);
    j = gd(2);
    if j > 1 && bf == 'l'
        j = j - 1;
    end
    if j < sm && bf == 'r'
        j = j + 1;
    end
    if i > 1 && bf == 'd'
        i = i - 1;
    end
    if i < st && bf == 'u'
        i = i + 1;
    end
    guidata(f,[i j]);
    N = Ms(j);
    set(f, 'Name', strcat("Eigenvector ", num2str(i), " Mesh ", num2str(j)));
    trisurf(N.TRIV, N.VERT(:, 1), N.VERT(:, 2), N.VERT(:, 3), N.evecs(:,track(i,j)));
    axis equal;
end

% Callback function for line plot
function Scroll_plot(~, ~, f, sm, track, bf)
    j = guidata(f);
    if j > 1 && bf == 'b'
        j = j - 1;
    end
    s = size(track); s = s(1);
    if j < s && bf == 'f'
        j = j + 1;
    end
    guidata(f,j);
    plot(1:sm, track(j,:));
    ylabel(strcat("Eigenfunction ", num2str(j)));
end