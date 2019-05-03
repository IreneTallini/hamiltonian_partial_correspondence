function [D, matches] = run_icp(M, N, r, C_init, max_iters)

C_init = C_init(:,1:r);

X = N.evecs(:,1:r)' * N.S;
Y = M.evecs' * M.S;
tree = ann('init', Y);

matches = ann('search', tree, C_init*X, 1)';

err = sum( sqrt(sum((C_init*X - Y(:,matches)).^2)) );
fprintf('(0) MSE: %f\n', err);

if max_iters == 0
    D = C_init;
    return
end

% Start iterations

[U,~,V] = svd(X * Y(:,matches)');
D = U * V(:,1:r)';
D = D';

matches = ann('search', tree, D*X, 1)';

err = sum( sqrt(sum((D*X - Y(:,matches)).^2)) );
fprintf('(1) MSE: %f\n', err);

ann('deinit', tree);
clear tree
