%input:
    %1) point_plane: punto sul piano
    %2) normal: normale al piano (coefficienti)
    %3) M: mesh da tagliare
%Output:
    %1) N: mesh tagliato
    %2) map_N_to_M: array in cui map_N_to_M(i) = vertice
    %   corrispondente in M
function [N, map_N_to_M] = cut_with_plane(point_plane, normal, M)
normal = normal/norm(normal);
is_over_plane = dot(M.VERT - point_plane, ...
    ones(M.n,1) * normal, 2) > 0;
N.VERT = M.VERT(is_over_plane, :);
N.n = sum(is_over_plane);
map_N_to_M = find(is_over_plane);
map_M_to_N = @(x) find(map_N_to_M == x);
is_over_plane_tri = all(ismember(M.TRIV, ...
    map_N_to_M), 2);
N.TRIV = M.TRIV(is_over_plane_tri, :);
N.m = sum(is_over_plane_tri);
N.TRIV = reshape(arrayfun(map_M_to_N, N.TRIV), ...
    N.m, 3);
end