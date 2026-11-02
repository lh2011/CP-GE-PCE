% Benchmark geometry from 
% "Inductance computation by incremental finite element analysis", 
% https://ieeexplore.ieee.org/document/767144
% rebuild in IGA

function [srf, regions, boundaries, params] = GEO_Transformer1(params)
    arguments
        params.preview (1,1) logical = true;
        params.a = 0.1;
        params.b = 0.1;
        params.c = 0.1;
        params.d = 0.1;
        params.e = 0.1;
        params.f = 0.1;
    end

    n_patches = 0;

    widths = [0, params.e, params.b, params.c, 2*params.d, params.c, params.b, params.e];
    cumWidths = cumsum(widths);
    heights = [0, params.e, params.a, params.e];
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
 
    % Patch definitions
    regions.Iron = [1:7,8,14,15,16:21];
    regions.IronInf = 11;
    regions.Phase1Out = 9;
    regions.Phase1In = 13;
    regions.Phase2Out = 10;
    regions.Phase2In = 12;

    boundaries.Dirichlet.patches = [1, 1, 2, 3, 4, 5, 6, 7, 7, 8, 14, 15, 15, 16, 17, 18, 19, 20, 21, 21];
    boundaries.Dirichlet.sides = [1, 3, 3, 3, 3, 3, 3, 2, 3, 1, 2, 1, 4, 4, 4, 4, 4, 4, 2, 4];

    if (params.preview)
        preview_nurbs(srf, regions);
        % for iPoint = 1:numel(point)
        %     scatter(point{iPoint}(1),point{iPoint}(2), "black", "filled");
        %     text(point{iPoint}(1)+0.001,point{iPoint}(2), string(iPoint))
        % end
        axis equal;
    end

    function addPatch(points)
        srf(end+1) = nrb4surf(point{points(1)}, point{points(2)}, point{points(3)}, point{points(4)});
        n_patches = n_patches + 1;
    end
end
