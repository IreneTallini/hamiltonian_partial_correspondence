function v = optimize_v(M,N,G,F,C,mu1,mu2,opt)
    % -------------------------------------------------------------------------
    % Solve for v
    % -------------------------------------------------------------------------
    A = C*N.evecs'*N.S*F;
    B = M.evecs'*M.S;

    areas = full(diag(M.S));
    target_area = full( sum(diag(N.S)) );


    x0 = M.evecs*((C*N.evecs')*(N.S*ones(N.n,1)));
  
    manifold = euclideanfactory(M.n,1);

    problem = {};
    problem.M = manifold;

    vfunc = mumford_shah(M.VERT, M.TRIV, M.S);
    problem.cost =  @(v) vfunc.cost(A, B, G, v, ones(size(v)), target_area, areas, mu1, mu2, opt);
    problem.egrad = @(v) vfunc.grad(A, B, G, v, ones(size(v)), target_area, areas, mu1, mu2, opt);

%     figure
%     checkgradient(problem);

    options.maxiter = 3e3;
    options.tolgradnorm = 1e-6;
    options.minstepsize = 1e-6;
    options.verbosity = 3;
    
    [v, cost, info, ~] = conjugategradient(problem, x0, options);
end
