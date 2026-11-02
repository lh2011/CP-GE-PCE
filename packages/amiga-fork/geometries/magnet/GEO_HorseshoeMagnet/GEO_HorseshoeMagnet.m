function [srf, regions, boundaries, params] = GEO_HorseshoeMagnet(params)
    arguments
        params.preview (1,1) logical = true;
        params.AirPadding = 2e-2;
        params.AirGapWidth = 2e-2;
        params.MagnetWidth = 5e-2;
        params.MagnetHeight = 2e-2;
        params.IronHeight = 5e-2;
        params.IronWidth = 2e-2;
    end
    
    n_patches = 0;
    widths = [0, params.AirPadding, params.IronWidth, params.MagnetWidth, params.AirGapWidth, params.IronWidth, params.AirPadding];
    cumWidths = cumsum(widths);
    heights = [0, params.AirPadding, params.MagnetHeight, params.IronHeight, params.MagnetHeight, params.AirPadding];
    cumHeights = cumsum(heights);
    
    for Hs = 1:numel(heights)
        for Ws = 1:numel(widths)
            point{Ws+(Hs-1)*numel(widths)} = [cumWidths(Ws), cumHeights(Hs)];
        end
    end
    
    temp = struct ('form', '', 'dim', 4, 'number', [], 'coefs', [], ...
        'knots', [], 'order', []);
    srf = temp([]);
    
    for Hs = 1:numel(heights)-1
        for Ws = 1:numel(widths)-1
            addPatch([Ws+(Hs-1)*numel(widths), Ws+1+(Hs-1)*numel(widths), ...
                Ws+(Hs)*numel(widths), Ws+1+(Hs)*numel(widths)])
        end
    end
    
    % region definitions
    regions.Air = [1:6, 7, 10, 12, 13, 15, 16, 18, 19, 22, 24, 25:30];
    regions.Iron = [8, 11, 14, 17, 20, 23];
    regions.MagnetLeft = 9;
    regions.MagnetRight = 21;
    % boundary definitions
    boundaries.Dirichlet.patches = [1,1,2,3,4,5,6,6,7,12,13,18,19,24,25,25,26,27,28,29,30,30];
    boundaries.Dirichlet.sides = [1,3,3,3,3,3,2,3,1,2,1,2,1,2,1,4,4,4,4,4,2,4];
    
    if (params.preview)
        preview_nurbs(srf, regions);
    
        % for i = 1:numel(point)
        %     scatter(point{i}(1), point{i}(2), "black", "filled");
        %     text(point{i}(1)+0.001, point{i}(2), string(i))
        % end
    end

    function addPatch(points)
        srf(end+1) = nrb4surf(point{points(1)}, point{points(2)}, point{points(3)}, point{points(4)});
        n_patches = n_patches + 1;
    end
end