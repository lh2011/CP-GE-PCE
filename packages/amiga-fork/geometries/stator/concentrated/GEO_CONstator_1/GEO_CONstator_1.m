% https://www.mathcha.io/editor/W05roux1SnLudOj1dGTKW8K0OcX2313ji9NEW5r

function [srf, regions, boundaries, params] = GEO_CONstator_1(params)
    arguments
        params.preview (1,1) logical = true;
        params.SLOTS = 6;
        params.Symmetry = 2;
    
        params.SRI = 47.8e-3/2;         % Stator Radius Inner
        params.SRO = 113e-3/2;          % Stator Radius Outer
        params.RAG = 23.7e-3;           % Radius Air Gap
        params.HY  = 11.6e-3;           % Height Yoke
        params.HSO = 0.98e-3;           % Height Slot Opening
        params.WT  = 14.8e-3;           % Width Tooth         
        params.WSO = 3.23e-3;           % Width Slot Opening
        params.LumpWindings (1,1) logical = false;
    end
    
    % calculated parameters
    beta = 2*pi/params.SLOTS;
    alpha = beta/2;
    gamma1 = asin(params.WSO/2 / params.SRI);
    
    L = cos(gamma1)*params.SRI + params.HSO;
    alpha2 = atan(cos(alpha)/(2*L/params.WT- sin(alpha)));
    alpha1 = alpha-alpha2;
    B = 2*L*tan(alpha1);
    
    delta = asin(params.WT/2/(params.SRO-params.HY));

    % define points
    point{1} = [0, 0];
    point{2} = params.SRI* [cos(alpha), -sin(alpha)];
    point{3} = params.SRI* [cos(gamma1), -sin(gamma1)];
    point{4} = point{3} + [params.HSO, 0];
    point{5} = point{4} + [0, params.WSO/2 - B/2];
    point{6} = (params.SRO - params.HY) * [cos(alpha-delta), -sin(alpha-delta)];
    point{7} = params.SRO * [cos(alpha), -sin(alpha)];
    gamma2 = atan2(-point{4}(2), point{4}(1));
    point{8} = (params.SRO - params.HY) * [cos(gamma2), -sin(gamma2)];
    point{9} = params.SRO * [cos(gamma2), -sin(gamma2)];
    gamma3 = atan2(-point{5}(2), point{5}(1));
    point{10} = params.SRI * [cos(gamma3), -sin(gamma3)];
    point{11} = [point{3}(1), 0];
    point{12} = [point{4}(1), 0];

    % acquire all angles
    angles = [];
    for iSrf=1:numel(point)
        angles(iSrf) = atan(point{iSrf}(2)/point{iSrf}(1));
    end

    % Air patches
    srf(1) = nrbruled(nrbcirc(params.RAG, point{1}, angles(2), angles(10)), nrbcirc(params.SRI, point{1}, angles(2), angles(10)));
    srf(2) = nrbruled(nrbcirc(params.RAG, point{1}, angles(10), angles(3)), nrbcirc(params.SRI, point{1}, angles(10), angles(3)));
    srf(3) = nrbruled(nrbcirc(params.RAG, point{1}, angles(3), 0), nrbcirc(params.SRI, point{1}, angles(3), 0));
    srf(4) = nrbruled(nrbcirc(params.SRI, point{1}, angles(3), 0), nrbline(point{4}, point{12}));
    
    % Iron parts
    srf(5) = nrbruled(nrbcirc(params.SRI, point{1}, angles(2), angles(10)), nrbcirc(vecnorm(point{5}), point{1}, angles(2), angles(10)));
    srf(6) = nrbruled(nrbcirc(params.SRI, point{1}, angles(10), angles(3)), nrbline(point{5}, point{4}));
    srf(7) = nrbruled(nrbcirc(vecnorm(point{5}), point{1}, angles(2), angles(10)), nrbline(point{7}, point{6}));
    srf(8) = nrbruled(nrbcirc(vecnorm(point{6}), point{1}, angles(6), angles(8)), nrbcirc(params.SRO, point{1}, angles(7), angles(9)));
    srf(9) = nrbruled(nrbcirc(vecnorm(point{8}), point{1}, angles(8), 0), nrbcirc(params.SRO, point{1}, angles(8), 0));
    
    % Windings
    srf(10) = nrbruled(nrbline(point{5}, point{4}), nrbcirc(vecnorm(point{8}), point{1}, angles(6), angles(9)));
    srf(11) = nrbruled(nrbline(point{4}, point{12}), nrbcirc(vecnorm(point{8}), point{1}, angles(8), 0));
    
    % Set same degree
    for iSrf = 1:numel(srf)
        srf(iSrf) = nrbdegelev(srf(iSrf), [3,3] - srf(iSrf).order);
    end

    % patch definitions
    regions.Iron = [5,6,7,8,9] ;
    regions.Air = [1,2,3,4];
    regions.Windings = [10,11];
    % Mirror half slot
    nsrf = numel(srf);
    mirror = diag(ones(4,1));
    mirror(2,2) = -1;
    for iSrf = 1:nsrf
        srf(iSrf+nsrf) = nrbtform(srf(iSrf), mirror);
    end
    regions.Iron = [regions.Iron, regions.Iron+nsrf];
    regions.Air = [regions.Air, regions.Air+nsrf];
    regions.Windings = [regions.Windings, regions.Windings+nsrf];
    regions.Slot1up = [21,22];
    regions.Slot1down = [10,11];
    
    % Rotate by half slot:
    nsrf = numel(srf);
    for iSrf = 1:nsrf
        srf(iSrf) = nrbtform(srf(iSrf), vecrotz(alpha));
    end

    boundaries.Dirichlet.patches = [8,9,19,20];
    boundaries.Dirichlet.sides = [4,4,4,4];
    boundaries.Coupling.patches = [1,2,3,12,13,14];
    boundaries.Coupling.sides = [3,3,3,3,3,3];

    % Copy and rotate
    nsrf = numel(srf);
    regionsInit = regions;
    boundariesInit = boundaries;
    for I = 1:params.SLOTS/params.Symmetry-1
        for iSrf = 1:nsrf
            srf(I*nsrf+iSrf) = nrbtform(srf(iSrf), vecrotz(I*2*pi/params.SLOTS));
        end
        regions.Iron = [regions.Iron, regionsInit.Iron+nsrf*I];
        regions.Air = [regions.Air, regionsInit.Air+nsrf*I];
        regions.("Slot"+num2str(I+1)+"up") = regionsInit.Slot1up + nsrf*I;
        regions.("Slot"+num2str(I+1)+"down") = regionsInit.Slot1down + nsrf*I;
        regions.Windings = [regions.Windings, regionsInit.Windings + nsrf*I];
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

    % Additional Boundaries if not full machine is simulated
    nsrf = numel(srf);
    if params.Symmetry ~= 1
        boundaries.Periodic1.patches = [1,5,7];
        boundaries.Periodic1.sides = [1,1,1];
        boundaries.Periodic2.patches = nsrf-[10,6,4];
        boundaries.Periodic2.sides = [1,1,1];
    end

    if (params.preview)
        preview_nurbs(srf, regions);
        % for iSrf = 1:numel(point)
        %     scatter(point{iSrf}(1), point{iSrf}(2), "black", "filled");
        %     text(point{iSrf}(1)+1e-4, point{iSrf}(2), string(iSrf))
        % end
    end
end