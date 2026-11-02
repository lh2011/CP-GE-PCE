% replica of 
% Friedrich, L. A. J. (2020). Isogeometric analysis of nonlinear eddy current problems.
% page 145
% https://doi.org/10.6100/5074ed73-fcec-4a93-8a5c-699abcf6e896
% value from "High-Order Methods Applied to Nonlinear Magnetostatic Problems" page 8

function [srf, regions, boundaries, params] = GEO_HorseshoeMagnet2(params)
    arguments
        params.preview (1,1) logical = true;
        params.Type {mustBeMember(params.Type, ["Magnet", "Coil"])} = "Magnet";
        params.HeightIronBar1 = 2e-2;
        params.HeightIronBar2 = 2e-2;
        params.HeightCoil = 2e-2;
        params.HeightAirGap = 1e-2;
        params.HeightMagnet = 3e-2;
        params.AirPadding = 5e-2;
        params.WidthIronBar = 2e-2;
        params.WidthMagnet = 5e-2;
    end
    
    n_patches = 0;
    widths = [0, params.AirPadding, params.WidthIronBar, params.WidthMagnet, params.WidthIronBar, params.AirPadding];
    cumWidths = cumsum(widths);
    heights = [0, params.AirPadding, params.HeightCoil, params.HeightMagnet, params.HeightCoil, params.HeightIronBar2, params.HeightAirGap, params.HeightIronBar1, params.AirPadding];
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
    if params.Type == "Magnet"
        regions.Iron = [12,14,17,19,22,24,32,33,34];
        regions.Air = [1:10, 11,15,16,18,20,21,23,25,26:30,31,35,36:40];
        regions.Magnet = 13;
    else
        regions.Iron = [12,13,14,17,19,22,24,32,33,34];
        regions.Air = [1:7,9,10, 11,15,16,20,21,23,25,26:30,31,35,36:40];
        regions.CoilIn = 8;
        regions.CoilOut = 18;
    end
    
    % boundary definitions
    boundaries.Dirichlet.patches = [1,1,2,3,4,5,5,6,10,11,15,16,20,21,25,26,30,31,35,36,36,37,38,39,40,40];
    boundaries.Dirichlet.sides = [1,3,3,3,3,2,3,1,2,1,2,1,2,1,2,1,2,1,2,1,4,4,4,4,2,4];

    if (params.preview)
        preview_nurbs(srf, regions);
    
        % for i = 1:numel(boundaries.Dirichlet.patches)
        %     nrbctrlplot(nrbextract(srf(boundaries.Dirichlet.patches(i)), boundaries.Dirichlet.sides(i)))
        % end
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