function [F1,F2] = calc_dense_descriptors(M, diam, options)

max_radius = diam * options.descr_max_radius/100;
vol_step = diam * options.intinv_volstep/100;
scales = linspace(vol_step, max_radius, 1+max_radius/vol_step);

[filterSizes, idx] = unique(round(scales./vol_step));
scales = scales(idx);
idx = logical(mod(filterSizes+1,2));
filterSizes = filterSizes(idx)+1;
scales = scales(idx);

if filterSizes(1)==1
    filterSizes = filterSizes(2:end);
    scales = scales(2:end);
end

F1 = calc_integral_invariant(M, vol_step, filterSizes, scales);
F2 = calc_shot(M.VERT', M.TRIV', 1:size(M.VERT,1), options.shot_bins, max_radius, options.shot_min_neighs)';

end
