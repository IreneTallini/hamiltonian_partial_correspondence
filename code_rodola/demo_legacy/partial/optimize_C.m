function [C, est_rank] = optimize_C(M, N, G, F, v, C_init, mu1, mu2)

k = size(M.evecs,2);

%% COMPUTE W MASK
est_rank = sum(N.evals-max(M.evals)<0);
W = zeros(k);
for i=1:k
    for j=1:k
        slope = est_rank/k;
        direction = [1 slope];
        direction = direction./norm(direction);
        W(i,j) = exp(-0.03*sqrt(i.^2 + j.^2))*norm(cross([direction 0], [i,j, 0]-[1 1 0]));        
    end
end
% imagesc(W);
d=zeros(1,k);
d(1:est_rank)=1;

if isempty(C_init)
    C_init = (max(max(W))-W)./max(max(W));
end

%% OPTIMIZE L2,1 (L2 at the end of file)

A = N.evecs'*N.S*F;
B = M.evecs'*M.S*diag(sparse(v))*G;

mu3 = 0;   % smoothness

manifold = euclideanfactory(k,k);
problem = {};

problem.M = manifold;

problem.cost = @(C) (...
	sum(sum((C*A-B).^2).^0.5) + ...
    mu1 * norm(C.*W,'fro')^2 + ...
    mu3 *(trace(C*diag(M.evals)*C') + trace(C'*diag(N.evals)*C)) + ...
    mu2 * (norm(C'*C,'fro')^2 - sum(diag(C'*C).^2) + sum((diag(C'*C) - d').^2) ));
%     mu2 * (norm(C'*C,'fro')^2 - sum((diag(C'*C).*d').^2)));

problem.egrad = @(C) (...
    norm_21_gradient(C,A,B) + ...
    mu1 * 2 * C.*W.*W + ...
    mu3 * 2 *(C*diag(M.evals) + diag(N.evals)*C) +...
    mu2 * 4*(C*C'*C - C.*repmat(sum(C.^2),k,1) + C.*(repmat(sum(C.^2),k,1) - repmat(d,k,1)) ));
%     mu2 * 4*(C*C'*C - C.*repmat(sum(C.^2),k,1).*repmat(d.^2,k,1)) );

%figure, checkgradient(problem)

options.maxiter = 1e4;
options.tolgradnorm = 1e-06;
options.minstepsize = 1e-06;
options.verbosity = 1;

%[C, xcost, info, options] = trustregions(problem);
[C, final_cost, ~, ~] = conjugategradient(problem,C_init,options);
fprintf('Final cost: %f\n', final_cost);

end
