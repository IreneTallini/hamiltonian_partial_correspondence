function desc = IntegralDesc(shape,volIndicator,volGrid,filtersizes,scales)

assert(length(scales)==length(filtersizes));

shape.X = shape.VERT(:,1);
shape.Y = shape.VERT(:,2);
shape.Z = shape.VERT(:,3);
n = length(shape.X);

% TODO - check smoothing filter size
dilateFilter = ones(5,5,5) / 5^3;
shapeG_dilated = double(convn(volIndicator,dilateFilter,'same'));

[shapeGrid.X,shapeGrid.Y,shapeGrid.Z] = meshgrid(volGrid{1},volGrid{2},volGrid{3}) ;

desc = zeros(n,numel(scales));

parfor i = 1:numel(scales)

    curRadius = scales(i);
    filterSize = filtersizes(i);
    
    assert(logical(mod(filterSize,2)));
    
    int = linspace(-curRadius, curRadius, filterSize);

    % create sphere filter
    [sphereX,sphereY,sphereZ] = meshgrid(int,int,int) ;
    sphereFilter = ((sphereX).^2 + (sphereY).^2 + (sphereZ).^2) <= (curRadius)^2;
    sphereFilter = sphereFilter / nnz(sphereFilter);

    % apply filter on volume matrix
    filterRes = convn(shapeG_dilated,sphereFilter,'same');
    filterRes = permute(filterRes,[2 1 3]);

    desc(:,i) = interp3(...
        shapeGrid.X, shapeGrid.Y, shapeGrid.Z, ...
        filterRes, ...
        shape.X, shape.Y, shape.Z, ...
        'linear');
    
%     drawme = filterRes;
%     figure
%     for j=1:size(drawme,3)
%         imagesc(drawme(:,:,j))
%         title(sprintf('Filter size: %d', filterSize))
%         axis equal
%         drawnow
%         pause(0.05)
%     end

end
