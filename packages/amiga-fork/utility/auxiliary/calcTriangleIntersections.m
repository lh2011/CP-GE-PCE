% function intersectionPoints = calcTriangleIntersections(triangles, nodes, limits)
function intersectionPoints = calcTriangleIntersections(v1, v2, v3, limits)
    % x/y | ntraingles | nlimits | trianglesides
    limits = reshape(limits, 1, []);
    intersectionPoints = NaN(2, size(v1, 1), numel(limits), 3);
    
    % v1 = nodes(triangles(:,1), :);
    % v2 = nodes(triangles(:,2), :);
    % v3 = nodes(triangles(:,3), :);
    
    lambda1 = (limits - v1(:, 3))./(v2(:,3) - v1(:,3));
    lambda1(lambda1<0) = NaN;
    lambda1(lambda1>1) = NaN;

    lambda2 = (limits - v2(:, 3))./(v3(:,3) - v2(:,3));
    lambda2(lambda2<0) = NaN;
    lambda2(lambda2>1) = NaN;

    lambda3 = (limits - v3(:, 3))./(v1(:,3) - v3(:,3));
    lambda3(lambda3<0) = NaN;
    lambda3(lambda3>1) = NaN;

    intersectionPoints(1, :, :, 1) = v1(:, 1) + (v2(:,1) - v1(:,1)).*lambda1;
    intersectionPoints(2, :, :, 1) = v1(:, 2) + (v2(:,2) - v1(:,2)).*lambda1;
    % pts(3, :, :, 1) = v1(:, 3) + (v2(:,3) - v1(:,3)).*lambda1;

    intersectionPoints(1, :, :, 2) = v2(:, 1) + (v3(:,1) - v2(:,1)).*lambda2;
    intersectionPoints(2, :, :, 2) = v2(:, 2) + (v3(:,2) - v2(:,2)).*lambda2;
    % pts(3, :, :, 2) = v2(:, 3) + (v3(:,3) - v2(:,3)).*lambda2;

    intersectionPoints(1, :, :, 3) = v3(:, 1) + (v1(:,1) - v3(:,1)).*lambda3;
    intersectionPoints(2, :, :, 3) = v3(:, 2) + (v1(:,2) - v3(:,2)).*lambda3;
    % pts(3, :, :, 3) = v3(:, 3) + (v1(:,3) - v3(:,3)).*lambda3;
end