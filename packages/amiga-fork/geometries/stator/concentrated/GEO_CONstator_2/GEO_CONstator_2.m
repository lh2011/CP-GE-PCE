

function [srf, regions, boundaries, params] = GEO_CONstator_2(params)
    arguments
        params.preview (1,1) logical = true;
        params.SLOTS = 12;
        params.Symmetry = 2;
    
        params.SRI = 20e-3;         % Stator Radius Inner
        params.SRO = 42e-3;         % Stator Radius Outer
        params.RAG = 42.5e-3;       % Radius Air Gap
        params.R1 = 1e-3;           % Radius 1
        params.R2 = 1.5e-3;         % Radius 2
        params.R3 = 1e-3;           % Radius 3
        params.HY = 4.5e-3;         % Height Yoke
        params.HSO = 1.6e-3;        % Height Slot Opening
        params.WT = 6.5e-3;         % Tooth Width
        params.WSO = 2.3e-3;        % Width Slot Opening
        params.LumpWindings (1,1) logical = false;
    end
    
    % calculated parameters
    center = [0,0,0];
    beta = pi/params.SLOTS;
    beta1 = atan(cos(beta)/((params.SRI+params.HY)/(0.5*params.WT) - sin(beta)));
    beta2 = beta - beta1;
    L = (params.SRI + params.HY)/cos(beta2);
    % define points
    point{1} = [cos(beta), sin(beta)] * (params.SRI + params.HY);
    point{2} = L*[cos(beta1), sin(beta1)];

    px = params.WSO/2/sin(beta); py = 0;
    a = 1; b = (2*px*cos(beta) + 2*py*sin(beta)); c = px^2 + py^2 - (params.SRO)^2;
    m = (-b + (b^2 - 4*a*c)^0.5)/(2*a);
    point{3} = [px, py] + (m- params.HSO)*[cos(beta), sin(beta)];

    point{4} = [point{3}(1), params.WT/2];
    [point{5}, point{6}, point{7}] = line2radius(point{1}, point{2}, point{4}, params.R3);

    [point{8}, point{9}, point{10}] = line2radius(point{6}, point{4}, point{3}, params.R2);
    point{11} = [params.SRI,0];
    point{12} = [(params.SRO + point{3}(1))/2, 0];
    % cross point line / arc
    px = point{3}(1); py = point{3}(2);
    a = 1; b = (2*px*cos(beta) + 2*py*sin(beta)); c = px^2 + py^2 - (params.SRO)^2;
    m = (-b + (b^2 - 4*a*c)^0.5)/(2*a);
    point{13} = point{3} + m*[cos(beta), sin(beta)];

    [point{14}, point{15}, point{16}] = line2radius(point{9}, point{3}, point{13}, params.R1);
    point{17} = [point{12}(1), point{16}(2)/2];
    point{18} = [cos(beta), sin(beta)] * vecnorm(point{7});
    
    point{19} = point{16} + vecnormalize(point{3} - point{16})*params.R1;
    point{20} = [cos(beta), sin(beta)]*(point{19}(1)*cos(beta) + point{19}(2)*sin(beta));
    point{21} = (point{7} + point{9} + point{20}) / 3;


    % surface definitions
    nrbRight = nrbglue(nrbdegelev(nrbline(point{1}, point{5}), 1), nrbcirc(params.R3, point{7}, pi+beta, 1.25*pi+beta/2));
    nrbRight.knots = nrbRight.knots/max(nrbRight.knots);
    nrbLeft = nrbreverse(nrbcirc(params.SRI, center, 0, beta));
    srf(1) = nrbruled(nrbLeft, nrbRight);
    srf(1) = nrbdegelev(srf(1), [0,1]);

    nrbTop1 = nrbcirc(params.R3, point{7}, 1.25*pi+beta/2, 1.5*pi);
    nrbTop2 = nrbline(point{6}, point{8});
    % nrbTop3 = nrbcirc(params.R2, point{10}, 1.5*pi, 1.75*pi);
    nrbTop3 = nrbcirc(params.R2, point{10}, 1.5*pi, 2*pi);
    nrbTop12 = nrbglue(nrbTop1, nrbdegelev(nrbTop2, 1));
    nrbTop = nrbglue(nrbTop12, nrbTop3);
    nrbTop.knots = nrbTop.knots/max(nrbTop.knots);
    nrbBot = nrbline(point{11}, point{12});
    srf(2) = nrbruled(nrbBot, nrbTop);
    srf(2) = nrbdegelev(srf(2), [0,1]);

    % nrbLeft1 = nrbcirc(params.R2, point{10}, 1.75*pi, 2*pi);
    nrbLeft2 = nrbdegelev(nrbline(point{9}, point{14}), 1);
    nrbLeft3 = nrbcirc(params.R1, point{16}, 0.75*pi+beta/2, pi);
    % nrbLeft12 = nrbglue(nrbLeft1, nrbLeft2);
    nrbLeft = nrbglue(nrbLeft2, nrbLeft3);
    nrbLeft.knots = nrbLeft.knots/max(nrbLeft.knots);
    srf(3) = nrbruled(nrbLeft, nrbdegelev(nrbline(point{12}, point{17}), 1));
    srf(3) = nrbdegelev(srf(3), [0,1]);

    nrbTop1 = nrbdegelev(nrbline(point{13}, point{15}), 1);
    nrbTop2 = nrbcirc(params.R1, point{16}, pi/2+beta, 0.75*pi+beta/2);
    nrbTop = nrbglue(nrbTop1, nrbTop2);
    nrbTop.knots = nrbTop.knots/max(nrbTop.knots);
    nrbTop = nrbreverse(nrbTop);
    nrbBot = nrbline(point{17}, point{17}/vecnorm(point{17})*params.SRO);
    nrbRight = nrbcirc(params.SRO, center, atan2(point{17}(2), point{17}(1)), atan2(point{13}(2), point{13}(1)));
    nrbLeft = nrbreverse(nrbextract(srf(3), 2));
    srf(4) = nrbcoons(nrbBot, nrbTop, nrbLeft, nrbRight);
    
    nrbTop = nrbextract(srf(4), 3);
    nrbLeft = nrbextract(srf(3), 4);
    nrbBot = nrbline(point{12}, [params.SRO, 0]);
    nrbRight = nrbcirc(params.SRO, center, 0, atan2(point{17}(2), point{17}(1)));
    srf(5) = nrbcoons(nrbBot, nrbTop, nrbLeft, nrbRight);

    srf(6) = nrbruled(nrbextract(srf(1), 4), nrbline(point{18}, point{7}));
    srf(6) = nrbdegelev(srf(6), [0,1]);
    srf(7) = nrbruled(nrbextract(srf(2), 4), nrbline(point{7}, point{21}));
    srf(7) = nrbdegelev(srf(7), [0,1]);
    srf(8) = nrbruled(nrbextract(srf(3), 3), nrbline(point{21}, point{20}));
    srf(8) = nrbdegelev(srf(8), [0,1]);

    nrbLeft = nrbreverse(nrbextract(srf(6), 4));
    nrbBot = nrbextract(srf(7), 4);
    nrbRight = nrbextract(srf(8), 4);
    nrbTop = nrbline(point{18}, point{20});
    srf(9) = nrbcoons(nrbBot, nrbTop, nrbLeft, nrbRight);

    nrbLeft = nrbextract(srf(8), 2);
    nrbBot = nrbextract(srf(4), 4);
    nrbRight = nrbcirc(params.SRO, center, atan2(point{13}(2), point{13}(1)), beta);
    nrbTop = nrbline(point{20}, [cos(beta), sin(beta)]*params.SRO);
    srf(10) = nrbcoons(nrbBot, nrbTop, nrbLeft, nrbRight);

    ang13 = atan2(point{13}(2), point{13}(1));
    ang17 = atan2(point{17}(2), point{17}(1));
    srf(11) = nrbruled(nrbextract(srf(5), 2), nrbcirc(params.RAG, center, 0, ang17));
    srf(12) = nrbruled(nrbextract(srf(4), 2), nrbcirc(params.RAG, center, ang17, ang13));
    srf(13) = nrbruled(nrbextract(srf(10), 2), nrbcirc(params.RAG, center, ang13, beta));

    % patch definitions
    regions.Iron = [1,2,3,4,5,14,15,16,17,18] ;
    regions.Windings = [6,7,8,9,19,20,21,22];
    regions.Air = [10,11,12,13,23,24,25,26];
    regions.Slot1up = [19,20,21,22];
    regions.Slot1down = [6,7,8,9];
    
    % Mirror half slot
    nsrf = numel(srf);
    mirror = diag(ones(4,1));
    mirror(2,2) = -1;
    for iSrf = 1:nsrf
        srf(iSrf+nsrf) = nrbtform(srf(iSrf), mirror);
    end
    % rotate by beta
    nsrf = numel(srf);
    for iSrf = 1:nsrf
        srf(iSrf) = nrbtform(srf(iSrf), vecrotz(beta));
    end

    boundaries.Dirichlet.patches = [1,14];
    boundaries.Dirichlet.sides = [3,3];
    boundaries.Coupling.patches = [11,12,13,24,25,26];
    boundaries.Coupling.sides = [4,4,4,4,4,4];

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
        boundaries.Periodic1.patches = [14,19,22,23,26];
        boundaries.Periodic1.sides = [1,1,4,4,2];
        boundaries.Periodic2.patches = nsrf-[25,20,17,16,13];
        boundaries.Periodic2.sides = [1,1,4,4,2];
    end

    if (params.preview)
        preview_nurbs(srf, regions);
        % for iSrf = 1:numel(point)
        %     scatter(point{iSrf}(1), point{iSrf}(2), "black", "filled");
        %     text(point{iSrf}(1)+1e-4, point{iSrf}(2), string(iSrf))
        % end
    end
end


function [P1, P2, Pr] = line2radius(p1, p2, p3, radius)

    sp1 = vecnormalize(p2-p1);
    sp2 = vecnormalize(p3-p2);
    
    alpha = acos(sum(sp1.*sp2));
    % If angle is too small, replace radius by line
    if alpha < 1e-4
        warning("Small radius angle! Radius is removed!");
    end
    d = radius*tan(alpha/2);
    P1 = p2-d*sp1;
    P2 = p2+d*sp2;
    rvec = vecnormalize(0.5*(sp1-sp2));
    Pr = p2-(radius^2+d^2)^0.5*rvec;
end