function D = calc_integral_invariant(M, volstep, filtersizes, scales)
    fprintf('Calculating integral invariants at %d dimensions... ', length(filtersizes));
    [G, ~, g, ~, ~] = mesh2volume_D(M, volstep);
    D = IntegralDesc(M, G, g, filtersizes, scales);
    fprintf('done.\n');
end
