function [evecs, evals] = dirichlet_efcts(M, inside, k)

assert(k<=length(inside))

[evecsR, evals] = eigs(M.W(inside, inside), M.A(inside, inside), k, 1e-10);
[evals, idx] = sort(diag(evals));

evecs = zeros(M.n, k);
evecs(inside, :) = evecsR(:, idx);

end

