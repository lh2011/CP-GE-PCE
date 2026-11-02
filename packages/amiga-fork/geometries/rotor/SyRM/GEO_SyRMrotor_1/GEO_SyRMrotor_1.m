
function [srf, regions, boundaries, params] = GEO_SyRMrotor_1(params)
    arguments
        params.preview (1,1) logical = true;
        params.POLES = 4;
        params.Symmetry = 4;
        
        params.RAG = 50.5e-3;                       % Radius Air Gap
        params.RRI = 15e-3;                         % Rotor Radius Inner
        params.RRO = 50e-3;                         % Rotor Radius Outer
        params.BN {mustBeGreaterThan(params.BN, 0), mustBeInteger} = 3;     % Barrier Number
        params.BW = [6e-3, 4e-3, 2e-3];             % Barrier Widths
        params.CW = [7e-3, 5e-3, 3e-3];             % Carrier Widths
        params.DRAD = 10e-3;                        % Distance first barrier from inner radius
        params.WOB = [2e-3, 1e-3, 0.5e-3];          % Width Outer Bridge
        params.WIB = [3e-3, 2e-3, 1e-3];            % Width Inner Bridge
        params.InnerBridge (1,1) logical = true;    % Generate inner bridges or set them to air?
    end
    % check if scalar values are given, in this case repeat the values
    if isscalar(params.BW)
        params.BW = repmat(params.BW, [1, params.BN]);
    end
    if isscalar(params.CW)
        params.CW = repmat(params.CW, [1, params.BN]);
    end
    if isscalar(params.WOB)
        params.WOB = repmat(params.WOB, [1, params.BN]);
    end
    if isscalar(params.WIB)
        params.WIB = repmat(params.WIB, [1, params.BN]);
    end
    % if barrier number is different from standard and not set: adjust standard parameters
    if numel(params.BW) ~= params.BN
        params.BW = ones(1, params.BN)*sin(pi/params.POLES)*params.RRO*0.3/params.BN;
    end
    if numel(params.CW) ~= params.BN
        params.CW = ones(1, params.BN)*sin(pi/params.POLES)*params.RRO*0.3/params.BN;
    end
    if numel(params.WOB) ~= params.BN
        params.WOB = ones(1, params.BN)*1e-3;
    end
    if numel(params.WIB) ~= params.BN
        params.WIB = ones(1, params.BN)*1e-3;
    end

    % calculated parameters
    beta = pi/params.POLES;
    % point definition
    for iBarrier = 1:params.BN
        bn = (iBarrier-1)*8;
        point{bn+1} = [params.RRI + params.DRAD + sum(params.CW(2:iBarrier)) + sum(params.BW(1:iBarrier-1)), 0];
        point{bn+2} = [point{bn+1}(1), params.WIB(iBarrier)/2];
        if iBarrier == 1
            point{bn+3} = [point{bn+2}(1), tan(beta)*point{bn+2}(1) - params.CW(1)/2/cos(beta)];
        else
            point{bn+3} = point{bn-1} + [1, tan(beta) - 1/cos(beta)]*params.CW(iBarrier);
        end

        px = point{bn+3}(1); py = point{bn+3}(2);
        a = 1; b = (2*px*cos(beta) + 2*py*sin(beta)); c = px^2 + py^2 - (params.RRO - params.WOB(iBarrier))^2;
        m = (-b + (b^2 - 4*a*c)^0.5)/(2*a);
        point{bn+4} = point{bn+3} + m*[cos(beta), sin(beta)];
        point{bn+5} = point{bn+1} + [params.BW(iBarrier), 0];
        point{bn+6} = point{bn+2} + [params.BW(iBarrier), 0];
        point{bn+7} = point{bn+3} + [1, tan(beta) - 1/cos(beta)]*params.BW(iBarrier);

        px = point{bn+7}(1); py = point{bn+7}(2);
        a = 1; b = (2*px*cos(beta) + 2*py*sin(beta)); c = px^2 + py^2 - (params.RRO - params.WOB(iBarrier))^2;
        m = (-b + (b^2 - 4*a*c)^0.5)/(2*a);
        point{bn+8} = point{bn+7} + m*[cos(beta), sin(beta)];
    end
    
    % region definitions, will be extended
    regions.Iron = [1,2,3];
    regions.Air = [4];
    % outer flux carrier
    center = [0,0,0];
    srf(1) = nrbruled(nrbcirc(params.RRI, center, beta/6, beta), nrbline(point{2}, point{3}));
    ang1 = atan2(point{4}(2), point{4}(1));
    ang2 = beta;
    srf(2) = nrbruled(nrbline(point{3}, [cos(beta), sin(beta)]*params.RRI), nrbcirc(params.RRO-params.WOB(1), center, ang1, ang2));        
    srf(3) = nrbruled(nrbcirc(params.RRO-params.WOB(1), center, ang1, ang2), nrbcirc(params.RRO, center, ang1, ang2));
    srf(4) = nrbruled(nrbcirc(params.RRO, center, ang1, ang2), nrbcirc(params.RAG, center, ang1, ang2));
    
    boundaries.Dirichlet.patches = [1];
    boundaries.Dirichlet.sides = [3];
    boundaries.Coupling.patches = [4];
    boundaries.Coupling.sides = [4];

    for iBarrier = 2:params.BN
        bn = (iBarrier-1)*8 - 4;
        % flux barrier 
        srf(bn+1) = nrbruled(nrbline(point{bn+2}, point{bn-2}), nrbline(point{bn+3}, point{bn-1}));
        ang1 = atan2(point{bn+4}(2), point{bn+4}(1));
        ang2 = atan2(point{bn}(2), point{bn}(1));
        srf(bn+2) = nrbruled(nrbline(point{bn+3}, point{bn-1}), nrbcirc(params.RRO-params.WOB(iBarrier-1), center, ang1, ang2));
        srf(bn+3) = nrbruled(nrbcirc(params.RRO-params.WOB(iBarrier-1), center, ang1, ang2), nrbcirc(params.RRO, center, ang1, ang2));
        % outer air
        srf(bn+4) = nrbruled(nrbcirc(params.RRO, center, ang1, ang2), nrbcirc(params.RAG, center, ang1, ang2));
        % 
        % % flux carrier 
        srf(bn+5) = nrbruled(nrbline(point{bn+6}, point{bn+2}), nrbline(point{bn+7}, point{bn+3}));
        ang1 = atan2(point{bn+8}(2), point{bn+8}(1));
        ang2 = atan2(point{bn+4}(2), point{bn+4}(1));
        srf(bn+6) = nrbruled(nrbline(point{bn+7}, point{bn+3}), nrbline(point{bn+8}, point{bn+4}));
        srf(bn+7) = nrbruled(nrbline(point{bn+8}, point{bn+4}), nrbcirc(params.RRO, center, ang1, ang2));
        % % outer air
        srf(bn+8) = nrbruled(nrbcirc(params.RRO, center, ang1, ang2), nrbcirc(params.RAG, center, ang1, ang2));
               
        regions.Air = [regions.Air, bn+1, bn+2, bn+4, bn+8];
        regions.Iron = [regions.Iron, bn+3, bn+5, bn+6, bn+7];
        boundaries.Coupling.patches = [boundaries.Coupling.patches, bn+4, bn+8];
        boundaries.Coupling.sides = [boundaries.Coupling.sides, 4, 4];
    end
    
    % last barrier
    bn = params.BN*8 - 4;
    % flux barrier 
    srf(bn+1) = nrbruled(nrbline(point{bn+2}, point{bn-2}), nrbline(point{bn+3}, point{bn-1}));
    ang1 = atan2(point{bn+4}(2), point{bn+4}(1));
    ang2 = atan2(point{bn}(2), point{bn}(1));
    srf(bn+2) = nrbruled(nrbline(point{bn+3}, point{bn-1}), nrbcirc(params.RRO-params.WOB(params.BN), center, ang1, ang2));
    srf(bn+3) = nrbruled(nrbcirc(params.RRO-params.WOB(params.BN), center, ang1, ang2), nrbcirc(params.RRO, center, ang1, ang2));
    % outer air
    srf(bn+4) = nrbruled(nrbcirc(params.RRO, center, ang1, ang2), nrbcirc(params.RAG, center, ang1, ang2));
    
    % last carrier
    ang1 = atan2(point{bn+2}(2), point{bn+2}(1));
    ang2 = atan2(point{bn+4}(2), point{bn+4}(1));
    srf(bn+5) = nrbruled(nrbline(point{bn+2}, point{bn+3}), nrbcirc(params.RRO-params.WOB(params.BN), center, ang1, ang2));
    srf(bn+6) = nrbruled(nrbcirc(params.RRO-params.WOB(params.BN), center, ang1, ang2), nrbcirc(params.RRO, center, ang1, ang2));
    srf(bn+7) = nrbruled(nrbcirc(params.RRO, center, ang1, ang2), nrbcirc(params.RAG, center, ang1, ang2));

    regions.Air = [regions.Air, bn+1, bn+2, bn+4, bn+7];
    regions.Iron = [regions.Iron, bn+3, bn+5, bn+6];
    boundaries.Coupling.patches = [boundaries.Coupling.patches, bn+4, bn+7];
    boundaries.Coupling.sides = [boundaries.Coupling.sides, 4, 4];

    % flux bridges
    nsrf = numel(srf);
    % first flux bridge
    srf(nsrf+1) = nrbruled(nrbcirc(params.RRI, center, 0, beta/6), nrbline(point{1}, point{2}));
    srf(nsrf+2) = nrbruled(nrbline(point{1}, point{2}), nrbline(point{5}, point{6}));
    regions.Iron = [regions.Iron, nsrf+1];
    if params.InnerBridge
        regions.Iron = [regions.Iron, nsrf+2];
    else
        regions.Air = [regions.Air, nsrf+2];
    end
    boundaries.Dirichlet.patches = [boundaries.Dirichlet.patches, nsrf+1];
    boundaries.Dirichlet.sides = [boundaries.Dirichlet.sides, 3];

    for iBarrier = 2:params.BN
        bn = (iBarrier-1)*8 - 4;
        srf(nsrf+iBarrier*2-1) = nrbruled(nrbline(point{bn+1}, point{bn+2}), nrbline(point{bn+5}, point{bn+6}));
        srf(nsrf+iBarrier*2) = nrbruled(nrbline(point{bn+5}, point{bn+6}), nrbline(point{bn+9}, point{bn+10}));
        regions.Iron = [regions.Iron, nsrf+iBarrier*2-1];
        if params.InnerBridge
            regions.Iron = [regions.Iron, nsrf+iBarrier*2];
        else
            regions.Air = [regions.Air, nsrf+iBarrier*2];
        end
    end
    bn = (params.BN-1)*8 - 4;
    % last iron piece
    nsrf = numel(srf);
    srf(nsrf+1) = nrbruled(nrbline(point{bn+9}, point{bn+10}), nrbcirc(params.RRO-params.WOB(params.BN), center, 0, ang1));
    srf(nsrf+2) = nrbruled(nrbcirc(params.RRO-params.WOB(params.BN), center, 0, ang1), nrbcirc(params.RRO, center, 0, ang1));
    srf(nsrf+3) = nrbruled(nrbcirc(params.RRO, center, 0, ang1), nrbcirc(params.RAG, center, 0, ang1));
    regions.Iron = [regions.Iron, nsrf+1, nsrf+2];
    regions.Air = [regions.Air, nsrf+3];
    boundaries.Coupling.patches = [boundaries.Coupling.patches, nsrf+3];
    boundaries.Coupling.sides = [boundaries.Coupling.sides, 4];
    
    % elevate to same degrees
    for iSrf = 1:numel(srf)
        srf(iSrf) = nrbdegelev(srf(iSrf), [3,3] - srf(iSrf).order);
    end

    % mirror patches
    nsrf = numel(srf);
    mirror = diag(ones(4, 1));
    mirror(2, 2) = -1;
    for iSrf = 1:nsrf
        srf(nsrf+iSrf) = nrbtform(srf(iSrf), mirror);
    end
    regions.Air = [regions.Air, regions.Air+nsrf];
    regions.Iron = [regions.Iron, regions.Iron+nsrf];

    boundaries.Dirichlet.patches = [boundaries.Dirichlet.patches, boundaries.Dirichlet.patches+nsrf];
    boundaries.Dirichlet.sides = [boundaries.Dirichlet.sides, boundaries.Dirichlet.sides];
    boundaries.Coupling.patches = [boundaries.Coupling.patches, boundaries.Coupling.patches+nsrf];
    boundaries.Coupling.sides = [boundaries.Coupling.sides, boundaries.Coupling.sides];
    if params.Symmetry ~= 1
        boundaries.Periodic1.patches = nsrf+[2,3,4];
        boundaries.Periodic1.sides = [2,2,2];
        nsrfHalf = nsrf;    % for second part of periodic BC
    end

    % rotate by half segment
    nsrf = numel(srf);
    for iSrf = 1:nsrf
        srf(iSrf) = nrbtform(srf(iSrf), vecrotz(beta));
    end

    % Copy and rotate
    nsrf = numel(srf);
    regionsInit = regions;
    boundariesInit = boundaries;
    for I = 1:params.POLES/params.Symmetry-1
        for i = 1:nsrf
            srf(I*nsrf+i) = nrbtform(srf(i), vecrotz(I*2*pi/params.POLES));
        end
        regions.Iron = [regions.Iron, regionsInit.Iron+nsrf*I];
        regions.Air = [regions.Air, regionsInit.Air+nsrf*I];
        boundaries.Dirichlet.patches = [boundaries.Dirichlet.patches, boundariesInit.Dirichlet.patches+nsrf*I];
        boundaries.Dirichlet.sides = [boundaries.Dirichlet.sides, boundariesInit.Dirichlet.sides];
        boundaries.Coupling.patches = [boundaries.Coupling.patches, boundariesInit.Coupling.patches+nsrf*I];
        boundaries.Coupling.sides = [boundaries.Coupling.sides, boundariesInit.Coupling.sides];
    end

    nsrf = numel(srf);
    if params.Symmetry ~= 1
        boundaries.Periodic2.patches = nsrf-2*nsrfHalf+[2,3,4];
        boundaries.Periodic2.sides = [2,2,2];
    end

    if (params.preview)
        preview_nurbs(srf, regions);
        % for i = 1:numel(point)
        %     scatter(point{i}(1), point{i}(2), "black", "filled");
        %     text(point{i}(1)+1e-4, point{i}(2), string(i))
        % end
        view(2)
        axis equal
    end
end