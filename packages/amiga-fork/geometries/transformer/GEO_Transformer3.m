function [srf, regions, boundaries, params] = GEO_Transformer3(params)
    arguments
        params.preview (1,1) logical  = true;
        params.w1 = 1e-2;
        params.w2 = 1e-2;
        params.h1 = 2e-2;
        params.h2 = 1e-2;
        params.WindingWidth = 0.4e-2;
        params.AirGap = 0.05e-2;
        params.WindingHeight = 1.8e-2;
    end
    
    regions.Iron = [];
    regions.Air = [];
    regions.Phase1in = [];
    regions.Phase1out = [];
    regions.Phase2in = [];
    regions.Phase2out = [];
    % Overwrite default parameter values

    AirGapTop = (params.h1 - params.WindingHeight) /2;

    n_patches = 0;

    % Points bottom 1
    point{1} = [0, 0];
    point{end+1} = [params.w1, 0];
    point{end+1} = [params.w1 + params.w2/3, 0];
    point{end+1} = [params.w1 + 2*params.w2/3, 0];
    point{end+1} = [params.w1 + params.w2, 0];
    point{end+1} = [2*params.w1 + params.w2, 0];
    % Points bottom 2
    point{end+1} = [0, params.h2];
    point{end+1} = [params.w1, params.h2];
    point{end+1} = [params.w1 + params.w2/3, params.h2];
    point{end+1} = [params.w1 + 2*params.w2/3, params.h2];
    point{end+1} = [params.w1 + params.w2, params.h2];
    point{end+1} = [2*params.w1 + params.w2, params.h2];
    % Points top 1
    point{end+1} = [0, params.h1+params.h2];
    point{end+1} = [params.w1, params.h1+params.h2];
    point{end+1} = [params.w1 + params.w2/3, params.h1+params.h2];
    point{end+1} = [params.w1 + 2*params.w2/3, params.h1+params.h2];
    point{end+1} = [params.w1 + params.w2, params.h1+params.h2];
    point{end+1} = [2*params.w1 + params.w2, params.h1+params.h2];
    % Points top 2
    point{end+1} = [0, params.h1+2*params.h2];
    point{end+1} = [params.w1, params.h1+2*params.h2];
    point{end+1} = [params.w1 + params.w2/3, params.h1+2*params.h2];
    point{end+1} = [params.w1 + 2*params.w2/3, params.h1+2*params.h2];
    point{end+1} = [params.w1 + params.w2, params.h1+2*params.h2];
    point{end+1} = [2*params.w1 + params.w2, params.h1+2*params.h2];
    % Points windings
    point{end+1} = [params.w1+params.AirGap, params.h2+AirGapTop];
    point{end+1} = [params.w1+params.AirGap, params.h1+params.h2-AirGapTop];
    point{end+1} = [params.w1+params.AirGap+params.WindingWidth, params.h2+AirGapTop];
    point{end+1} = [params.w1+params.AirGap+params.WindingWidth, params.h1+params.h2-AirGapTop];
    point{end+1} = [params.w1+params.w2-params.AirGap, params.h2+AirGapTop];
    point{end+1} = [params.w1+params.w2-params.AirGap, params.h1+params.h2-AirGapTop];
    point{end+1} = [params.w1+params.w2-params.AirGap-params.WindingWidth, params.h2+AirGapTop];
    point{end+1} = [params.w1+params.w2-params.AirGap-params.WindingWidth, params.h1+params.h2-AirGapTop];
    
    % Points windings outer left
    point{end+1} = [-params.AirGap, params.h2+AirGapTop];
    point{end+1} = [-params.AirGap, params.h1+params.h2-AirGapTop];
    point{end+1} = [-params.WindingWidth-params.AirGap, 0];
    point{end+1} = [-params.WindingWidth-params.AirGap, params.h2];
    point{end+1} = [-params.WindingWidth-params.AirGap, params.h2+AirGapTop];
    point{end+1} = [-params.WindingWidth-params.AirGap, params.h1+params.h2-AirGapTop];
    point{end+1} = [-params.WindingWidth-params.AirGap, params.h1+params.h2];
    point{end+1} = [-params.WindingWidth-params.AirGap, params.h1+2*params.h2];
    % Points windings outer right
    point{end+1} = [2*params.w1+1*params.w2+params.AirGap, params.h2+AirGapTop];
    point{end+1} = [2*params.w1+1*params.w2+params.AirGap, params.h1+params.h2-AirGapTop];
    point{end+1} = [2*params.w1+1*params.w2+params.WindingWidth+params.AirGap, 0];
    point{end+1} = [2*params.w1+1*params.w2+params.WindingWidth+params.AirGap, params.h2];
    point{end+1} = [2*params.w1+1*params.w2+params.WindingWidth+params.AirGap, params.h2+AirGapTop];
    point{end+1} = [2*params.w1+1*params.w2+params.WindingWidth+params.AirGap, params.h1+params.h2-AirGapTop];
    point{end+1} = [2*params.w1+1*params.w2+params.WindingWidth+params.AirGap, params.h1+params.h2];
    point{end+1} = [2*params.w1+1*params.w2+params.WindingWidth+params.AirGap, params.h1+2*params.h2];

    % Bottom iron
    srf(1) = nrb4surf(point{1}, point{2}, point{7}, point{8});
    n_patches = n_patches + 1;
    regions.Iron = [regions.Iron, n_patches];    
    for p = 1:4
        addPatch([1+p,2+p,7+p,8+p], 'Iron')
    end
    % Center iron
    addPatch([7,8,13,14], 'Iron')
    addPatch([11,12,17,18], 'Iron')
    for p = 1:5
        addPatch([12+p,13+p,18+p,19+p], 'Iron')
    end
    % AIR left trafo
    addPatch([35,1,36,7], 'Air')
    addPatch([36,7,37,33], 'Air')
    addPatch([38,34,39,13], 'Air')
    addPatch([39,13,40,19], 'Air')
    addPatch([33,7,34,13], 'Air')

    % Air center trafo
    addPatch([8,9,25,27], 'Air')
    addPatch([26,28,14,15], 'Air')
    addPatch([9,10,27,31], 'Air')
    addPatch([27,31,28,32], 'Air')
    addPatch([28,32,15,16], 'Air')
    addPatch([10,11,31,29], 'Air')
    addPatch([32,30,16,17], 'Air')
    addPatch([8,25,14,26], 'Air')
    addPatch([29,11,30,17], 'Air')

    % Air right trafo
    addPatch([6,43,12,44], 'Air')
    addPatch([12,44,41,45], 'Air')
    addPatch([42,46,18,47], 'Air')
    addPatch([18,47,24,48], 'Air')
    addPatch([12,41,18,42], 'Air')
     
    % Windings
    addPatch([37,33,38,34], 'Phase1in')
    addPatch([25,27,26,28], 'Phase1out')
    addPatch([31,29,32,30], 'Phase2in')
    addPatch([41,45,42,46], 'Phase2out')

    boundaries.Dirichlet.patches = [1,2,3,4,5,8,9,10,11,12,13,13,14,15,16,16,27,27,28,29,30,30,32,35];
    boundaries.Dirichlet.sides = [3,3,3,3,3,4,4,4,4,4,1,3,1,1,1,4,2,3,2,2,2,4,1,2];
    
    if (params.preview)
        preview_nurbs(srf, regions);
        for iPoint = 1:numel(point)
            scatter(point{iPoint}(1),point{iPoint}(2), "black", "filled");
            text(point{iPoint}(1)+0.001,point{iPoint}(2), string(iPoint))
        end
    end

    function addPatch(points, material)
        srf(end+1) = nrb4surf(point{points(1)},point{points(2)},point{points(3)},point{points(4)});
        n_patches = n_patches + 1;
        regions.(material) = [regions.(material), n_patches];
    end
end