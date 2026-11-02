% TBD is this COMPUMAG TEAM benchmark ?
function [srf, regions, boundaries, params] = GEO_Transformer2(params)
    arguments
        params.preview (1,1) logical = true;
        params.L1 = 30e-3;
        params.L2 = 21e-3;
        params.L3 = 63e-3;
        params.L4 = 54e-3;
        params.H1 = 105e-3;
        params.H2 = 120e-3;
        params.H3 = 150e-3;
        params.H4 = 225e-3;
    end
    
    L1 = params.L1; L2 = params.L2; L3 = params.L3; L4 = params.L4;
    H1 = params.H1; H2 = params.H2; H3 = params.H3; H4 = params.H4;

    n_patches = 0;

    widths = [0,L4,L2,L1,L2,L3,L2,L1,L2,L3,L2,L1,L2,L4];
    cumWidths = cumsum(widths);
    cumWidths = cumWidths - cumWidths(end)/2; % center
    heights = [0,H1,H2-H1,H3-H2,H4-H3];
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
    regions.Iron = [3,7,11,16,20,24,29:37];
    regions.Air = [1,5,9,13,14,15,17:19,21:23,25,26,27,28,38,39,40:52];
    regions.Phase1up = [2];
    regions.Phase1down = [4];
    regions.Phase2up = [6];
    regions.Phase2down = [8];
    regions.Phase3up = [10];
    regions.Phase3down = [12];

    regions.Windings = [2,4,6,8,10,12];

    boundaries.Dirichlet.patches = [1,13,14,26,27,39,40,40:52,52];
    boundaries.Dirichlet.sides = [1,2,1,2,1,2,1,repmat(4,1,13),2];
    
    if (params.preview)
        preview_nurbs(srf, regions)
        % for iPoint = 1:numel(point)
        %     scatter(point{iPoint}(1),point{iPoint}(2), "black", "filled");
        %     text(point{iPoint}(1)+0.001,point{iPoint}(2), string(iPoint))
        % end
    end

    function addPatch(points)
        srf(end+1) = nrb4surf(point{points(1)}, point{points(2)}, point{points(3)}, point{points(4)});
        n_patches = n_patches + 1;
    end
end
