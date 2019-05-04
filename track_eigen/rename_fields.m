function N = rename_fields(M, dataset_name)
if (dataset_name == "SGP")
    N.VERT = M.xyz;
    N.TRIV = M.tri;
    n = size(N.VERT);
    N.n = n(1);
    m = size(N.TRIV);
    N.m = m(1);
    N.gt = M.gt;
end

if (dataset_name == "TOSCA")
    N.VERT(:, 1) = M.X;
    N.VERT(:, 2) = M.Y;
    N.VERT(:, 3) = M.Z;
    N.TRIV = M.TRIV;
    n = size(N.VERT);
    N.n = n(1);
    m = size(N.TRIV);
    N.m = m(1);
end
end