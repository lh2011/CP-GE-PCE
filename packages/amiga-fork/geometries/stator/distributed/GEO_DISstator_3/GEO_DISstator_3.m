% 

function [srf, regions, boundaries, params] = GEO_DISstator_3(params)
    arguments
        params.preview = true;
        params.Symmetry = 8;
        params.SLOTS = 48;              % Nuber of stator slots
        params.RAG = 78.9e-3;           % Radius Air Gap
        params.SRI = 79.1e-3;             % Stator Radius Inner
        params.SRO = 116e-3;            % Stator Radius Outer
        
        params.HSO = 1.1e-3;            % Height Slot Opening
        params.HY = 17.1e-3;            % Height Yoke
        params.WSO = 2.5e-3;            % Width Slot Opening
        params.WS = 5.1e-3;             % Width Slot
        params.LumpWindings (1,1) logical = false;
    end
    
    % patch definitions
    regions.Iron = [4,9,10,11,12,16,21,22,23,24];
    regions.Air = [1,2,3,13,14,15];
    regions.Slot1up = [7,8,19,20];
    regions.Slot1down = [5,6,17,18];
    regions.Windings = [regions.Slot1up, regions.Slot1down];
    
    % calculated parameters
    a1 = asin((params.WSO/2)/(params.SRI));
    a2 = pi/params.SLOTS;
    % a2 = asin((params.WSO/2)/(params.SRI+params.HSO));

    point{1} =  [cos(a1)*params.SRI+params.HSO, 0];
    point{2} =  [cos(a1)*params.SRI+params.HSO, params.WSO/2];
    point{3} =  [cos(a1)*params.SRI+params.HSO, params.WS/2];
    LS = params.SRO-params.HY-point{1}(1);
    point{4} =  point{1} + [LS/2, 0];
    point{5} =  point{2} + [LS/2, 0];
    point{6} =  point{3} + [LS/2, 0];
    point{7} =  point{1} + [LS, 0];
    point{8} =  point{2} + [LS, 0];
    point{9} =  point{3} + [LS, 0];

    % Air patches
    srf(1) = nrbruled(nrbcirc(params.RAG,[0,0],0,a1),nrbcirc(params.SRI,[0,0],0,a1));
    srf(2) = nrbruled(nrbcirc(params.RAG,[0,0],a1,a2),nrbcirc(params.SRI,[0,0],a1,a2));
    srf(3) = nrbruled(nrbcirc(params.SRI,[0,0],0,a1), nrbline(point{1}, point{2}));
    % Iron patch
    srf(4) = nrbruled(nrbcirc(params.SRI,[0,0],a1,a2), nrbline(point{2}, point{3}));
    % Coil patches
    srf(5) = nrbruled(nrbextract(srf(3), 4), nrbline(point{4}, point{5}));
    srf(6) = nrbruled(nrbextract(srf(4), 4), nrbline(point{5}, point{6}));
    srf(7) = nrbruled(nrbextract(srf(5), 4), nrbline(point{7}, point{8}));
    srf(8) = nrbruled(nrbextract(srf(6), 4), nrbline(point{8}, point{9}));
    % Iron patches
    srf(9) = nrbruled(nrbreverse(nrbextract(srf(4),2)), nrbline(point{6}, [cos(a2), sin(a2)]*vecnorm(point{6})));
    srf(10) = nrbruled(nrbextract(srf(9),4), nrbline(point{9}, [cos(a2), sin(a2)]*params.SRO));
    srf(11) = nrbruled(nrbextract(srf(7),4), nrbcirc(params.SRO, [0,0], 0, a1));
    srf(12) = nrbruled(nrbextract(srf(8),4), nrbcirc(params.SRO, [0,0], a1, a2));
    
    % mirror first half winding
    nsrf = numel(srf);
    mirror = diag(ones(4,1));
    mirror(2,2) = -1;
    for iSrf = 1:nsrf
        srf(nsrf+iSrf) = nrbtform(srf(iSrf), mirror);
    end
    % Rotate again by half a slot and upward facing
    nsrf = numel(srf);
    for iSrf = 1:nsrf
        srf(iSrf) = nrbtform(srf(iSrf), vecrotz(pi/params.SLOTS));
    end

    % Boundaries
    boundaries.Dirichlet.patches = [11,12,23,24];
    boundaries.Dirichlet.sides = [4,4,4,4];

    boundaries.Coupling.patches = [1,2,13,14];
    boundaries.Coupling.sides = [3,3,3,3];
    
    % Copy and rotate
    nsrf = numel(srf);
    IronInit = regions.Iron;
    AirInit = regions.Air;
    WindingsInit = regions.Windings;
    boundariesInit = boundaries;
    for I = 1:params.SLOTS/params.Symmetry-1
        for iSrf = 1:nsrf
            srf(I*nsrf+iSrf) = nrbtform(srf(iSrf), vecrotz(I*2*pi/params.SLOTS));
        end
        regions.Iron = [regions.Iron, IronInit+nsrf*I];
        regions.Air = [regions.Air, AirInit+nsrf*I];
        regions.Windings = [regions.Windings, WindingsInit+nsrf*I];
        regions.("Slot"+num2str(I+1)+"up") = regions.Slot1up + nsrf*I;
        regions.("Slot"+num2str(I+1)+"down") = regions.Slot1down + nsrf*I;
        boundaries.Dirichlet.patches = [boundaries.Dirichlet.patches, boundariesInit.Dirichlet.patches+nsrf*I];
        boundaries.Dirichlet.sides = [boundaries.Dirichlet.sides, boundariesInit.Dirichlet.sides];
        boundaries.Coupling.patches = [boundaries.Coupling.patches, boundariesInit.Coupling.patches+nsrf*I];
        boundaries.Coupling.sides = [boundaries.Coupling.sides, boundariesInit.Coupling.sides];
    end
    if params.LumpWindings
        regionsNoSlots.Iron = regions.Iron;
        regionsNoSlots.Air = regions.Air;
        regionsNoSlots.Windings = regions.Windings;
        regions = regionsNoSlots;
    else
        regions = rmfield(regions, "Windings");
    end

    nsrf = numel(srf);
    % Additional Boundaries if not full machine is simulated
    if params.Symmetry ~= 1
        boundaries.Periodic1.patches = [14,21,22];
        boundaries.Periodic1.sides = [2,2,2];
        boundaries.Periodic2.patches = nsrf-[22,15,14];
        boundaries.Periodic2.sides = [2,2,2];
    end

    if (params.preview)
        preview_nurbs(srf, regions);
    end
end
