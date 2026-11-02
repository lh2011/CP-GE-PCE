
function [srf, regions, boundaries, params] = GEO_SyRMrotor_2(params)
    arguments
        params.preview (1, 1) logical = true;
        params.POLES = 4;
        params.Symmetry = 4;
        
        params.RAG = 50.5e-3;                       % Radius Air Gap
        params.RRI = 25e-3;                         % Rotor Radius Inner
        params.RRO = 50e-3;                         % Rotor Radius Outer
        params.BN {mustBeGreaterThan(params.BN, 0), mustBeInteger} = 3;     % Barrier Number
        params.AngleBarrier (1,:) double = [];   % Barrier Angles (array with values for each barrier)
        params.AngleCarrier (1,:) double = [];   % Carrier Angles (array with values for each carrier)
        params.CW (1,:) double = [];             % Carrier Width at center (array with values for each barrier)
        params.BW (1,:) double = [];             % Barrier Width at center (array with values for each carrier)
        params.WOB (1,:) double = [];            % Width Outer Bridge (array with values for each carrier)
        params.R1 (1,:) double = [];             % Radius 1 (array with values for each carrier)
        params.R2 (1,:) double = [];             % Radius 2 (array with values for each carrier)
    end

    if isscalar(params.AngleBarrier)
        params.AngleBarrier = repmat(params.AngleBarrier, [1, params.BN]);
    end
    if isscalar(params.AngleCarrier)
        params.AngleCarrier = repmat(params.AngleCarrier, [1, params.BN]);
    end
    if isscalar(params.CW)
        params.CW = repmat(params.CW, [1, params.BN]);
    end
    if isscalar(params.BW)
        params.BW = repmat(params.BW, [1, params.BN]);
    end
    if isscalar(params.WOB)
        params.WOB = repmat(params.WOB, [1, params.BN]);
    end
    if isscalar(params.R1)
        params.R1 = repmat(params.R1, [1, params.BN]);
    end
    if isscalar(params.R2)
        params.R2 = repmat(params.R2, [1, params.BN]);
    end

    % calculated parameters
    beta = pi/params.POLES;
    % if barrier number is different from standard and not set: adjust standard parameters
    if numel(params.AngleBarrier) ~= params.BN
        params.AngleBarrier = ones(1, params.BN)*beta/(params.BN+1)/2.3;
    end
    if numel(params.AngleCarrier) ~= params.BN
        params.AngleCarrier = ones(1, params.BN)*beta/(params.BN+1)/2.3;
    end
    if numel(params.WOB) ~= params.BN
        params.WOB = ones(1, params.BN)*5e-4;
    end
    if numel(params.CW) ~= params.BN
        params.CW = ones(1, params.BN)*(params.RRO-params.RRI-params.WOB(end))/(params.BN+0.5)/2;
    end
    if numel(params.BW) ~= params.BN
        params.BW = ones(1, params.BN)*(params.RRO-params.RRI-params.WOB(end))/(params.BN+0.5)/2;
    end
    if numel(params.R1) ~= params.BN
        params.R1 = ones(1, params.BN)*1e-3/params.BN;
    end
    if numel(params.R2) ~= params.BN
        params.R2 = ones(1, params.BN)*1e-3/params.BN;
    end


    point = cell(10, params.BN);
    % point definition
    for iBarrier = 1:params.BN
        % circle 1
        point{1, iBarrier} = [params.RRI + sum(params.CW(1:iBarrier)) + sum(params.BW(1:(iBarrier-1))), 0];
        angle = beta - sum(params.AngleBarrier(1:iBarrier)) - sum(params.AngleCarrier(1:(iBarrier-1)));
        point{2, iBarrier} = [cos(angle), sin(angle)]*(params.RRO-params.WOB(iBarrier));
        
        alpha = atan2(point{2, iBarrier}(2), point{2, iBarrier}(1)-point{1, iBarrier}(1));
        L = vecnorm(point{2, iBarrier}-point{1, iBarrier});
        R1 = L/(2*cos(alpha));
        point{3, iBarrier} = point{1, iBarrier} + [R1, 0];    % center circle 1

        % edge fillet 1
        d1 = vecnorm(point{3, iBarrier});
        rho1 = params.RRO - params.WOB(iBarrier) - params.R1(iBarrier);
        rho2 = R1 - params.R1(iBarrier);
        a = (rho1^2 - rho2^2 + d1 ^2)/(2*d1);
        h = (rho1^2- a^2)^0.5;
        u = point{3, iBarrier} / vecnorm(point{3, iBarrier});
        n = [-u(2), u(1)];
        
        point{4, iBarrier} = a*u + h*n;

        point{5, iBarrier} = point{4, iBarrier}/vecnorm(point{4, iBarrier}) * (params.RRO - params.WOB(iBarrier));
        point{6, iBarrier} = point{3, iBarrier} + (point{4, iBarrier}-point{3, iBarrier})/vecnorm(point{4, iBarrier}-point{3, iBarrier}) * R1;
        n = (point{6, iBarrier}-point{4, iBarrier}) + (point{5, iBarrier}-point{4, iBarrier});
        n = n/vecnorm(n);
        point{7, iBarrier} = point{4, iBarrier} + params.R1(iBarrier)*n;

        point{8, iBarrier} = point{1, iBarrier} + [1/3*params.BW(iBarrier), 0];

        alpha = atan2(point{4, iBarrier}(2), point{4, iBarrier}(1)-point{8, iBarrier}(1));
        L = vecnorm(point{4, iBarrier}-point{8, iBarrier});
        R1 = L/(2*cos(alpha));
        point{9, iBarrier} = point{8, iBarrier} + [R1, 0];

        a1shift = atan2(point{7, iBarrier}(2)-point{3, iBarrier}(2), point{7, iBarrier}(1)-point{3, iBarrier}(1));
        a2shift = atan2(point{6, iBarrier}(2)-point{3, iBarrier}(2), point{6, iBarrier}(1)-point{3, iBarrier}(1));
        ashift = a2shift - a1shift;
        a1 = atan2(point{4, iBarrier}(2)-point{9, iBarrier}(2), point{4, iBarrier}(1)-point{9, iBarrier}(1));
        
        point{10, iBarrier} = point{9, iBarrier} + [cos(a1+ashift), sin(a1+ashift)] * vecnorm(point{4, iBarrier} - point{9, iBarrier});
        point{11, iBarrier} = point{9, iBarrier} + [cos((a1+ashift+pi)/2), sin((a1+ashift+pi)/2)] * vecnorm(point{4, iBarrier} - point{9, iBarrier});

        point{12, iBarrier} = point{2, iBarrier}/vecnorm(point{2, iBarrier})*params.RRO;
        point{13, iBarrier} = point{2, iBarrier}/vecnorm(point{2, iBarrier})*params.RAG;

        a1 = atan2(point{6, iBarrier}(2)-point{3, iBarrier}(2), point{6, iBarrier}(1)-point{3, iBarrier}(1));
        point{14, iBarrier} = point{3, iBarrier} + [cos((a1+pi)/2), sin((a1+pi)/2)] * vecnorm(point{6, iBarrier} - point{3, iBarrier});


        % circle 2 %%%%%%%%%%%%%%%%%%%%
        point{21, iBarrier} = [params.RRI + sum(params.CW(1:iBarrier)) + sum(params.BW(1:iBarrier)), 0];
        angle = beta - sum(params.AngleBarrier(1:iBarrier)) - sum(params.AngleCarrier(1:iBarrier));
        point{22, iBarrier} = [cos(angle), sin(angle)]*(params.RRO-params.WOB(iBarrier));
        
        alpha = atan2(point{22, iBarrier}(2), point{22, iBarrier}(1)-point{21, iBarrier}(1));
        L = vecnorm(point{22, iBarrier}-point{21, iBarrier});
        R1 = L/(2*cos(alpha));
        point{23, iBarrier} = point{21, iBarrier} + [R1, 0];    % center circle 1

        % edge fillet 2
        d1 = vecnorm(point{23, iBarrier});
        rho1 = params.RRO - params.WOB(iBarrier) - params.R2(iBarrier);
        rho2 = R1 + params.R2(iBarrier);
        a = (rho1^2 - rho2^2 + d1 ^2)/(2*d1);
        h = (rho1^2- a^2)^0.5;
        u = point{23, iBarrier} / vecnorm(point{23, iBarrier});
        n = [-u(2), u(1)];
        
        point{24, iBarrier} = a*u + h*n;

        point{25, iBarrier} = point{24, iBarrier}/vecnorm(point{24, iBarrier}) * (params.RRO - params.WOB(iBarrier));
        point{26, iBarrier} = point{23, iBarrier} + (point{24, iBarrier}-point{23, iBarrier})/vecnorm(point{24, iBarrier}-point{23, iBarrier}) * R1;
        n = (point{26, iBarrier}-point{24, iBarrier}) + (point{25, iBarrier}-point{24, iBarrier});
        n = n/vecnorm(n);
        point{27, iBarrier} = point{24, iBarrier} + params.R2(iBarrier)*n;

        point{28, iBarrier} = point{1, iBarrier} + [2/3*params.BW(iBarrier), 0];

        alpha = atan2(point{24, iBarrier}(2), point{24, iBarrier}(1)-point{28, iBarrier}(1));
        L = vecnorm(point{24, iBarrier}-point{28, iBarrier});
        R1 = L/(2*cos(alpha));
        point{29, iBarrier} = point{28, iBarrier} + [R1, 0];

        a1shift = atan2(point{27, iBarrier}(2)-point{23, iBarrier}(2), point{27, iBarrier}(1)-point{23, iBarrier}(1));
        a2shift = atan2(point{26, iBarrier}(2)-point{23, iBarrier}(2), point{26, iBarrier}(1)-point{23, iBarrier}(1));
        ashift = a2shift - a1shift;
        a1 = atan2(point{24, iBarrier}(2)-point{29, iBarrier}(2), point{24, iBarrier}(1)-point{29, iBarrier}(1));
        
        point{30, iBarrier} = point{29, iBarrier} + [cos(a1+ashift), sin(a1+ashift)] * vecnorm(point{24, iBarrier} - point{29, iBarrier});
        
        point{31, iBarrier} = point{29, iBarrier} + [cos((a1+ashift+pi)/2), sin((a1+ashift+pi)/2)] * vecnorm(point{24, iBarrier} - point{29, iBarrier});

        point{32, iBarrier} = point{22, iBarrier}/vecnorm(point{22, iBarrier})*params.RRO;
        point{33, iBarrier} = point{22, iBarrier}/vecnorm(point{22, iBarrier})*params.RAG;

        a1 = atan2(point{26, iBarrier}(2)-point{23, iBarrier}(2), point{26, iBarrier}(1)-point{23, iBarrier}(1));
        point{34, iBarrier} = point{23, iBarrier} + [cos((a1+pi)/2), sin((a1+pi)/2)] * vecnorm(point{26, iBarrier} - point{23, iBarrier});

    end

    regions.Iron = [];
    regions.Air = [];
    boundaries.Dirichlet.patches = [];
    boundaries.Dirichlet.sides = [];
    boundaries.Coupling.patches = [];
    boundaries.Coupling.sides = [];

    for iBarrier = 1:params.BN
        % flux barrier 1
        nrb1 = threePointCircle(point{7, iBarrier}, point{6, iBarrier}, point{4, iBarrier});
        nrb2 = threePointCircle(point{6, iBarrier}, point{14, iBarrier}, point{3, iBarrier});
        nrbLeft = nrbglue(nrb1, nrb2, 2, 1);
        nrbLeft.knots = nrbLeft.knots / max(nrbLeft.knots);

        nrb1 = threePointCircle(point{4, iBarrier}, point{10, iBarrier}, point{9, iBarrier});
        nrb2 = threePointCircle(point{10, iBarrier}, point{11, iBarrier}, point{9, iBarrier});
        nrbRight = nrbglue(nrb1, nrb2, 2, 1);
        nrbRight.knots = nrbRight.knots / max(nrbRight.knots);

        srf(1, iBarrier) = nrbruled(nrbLeft, nrbRight);

        % flux barrier 2
        nrb1 = threePointCircle(point{26, iBarrier}, point{27, iBarrier}, point{24, iBarrier});
        nrb2 = threePointCircle(point{26, iBarrier}, point{34, iBarrier}, point{23, iBarrier});
        nrbLeft = nrbglue(nrb1, nrb2, 1, 1);
        nrbLeft.knots = nrbLeft.knots / max(nrbLeft.knots);
        nrbLeft = nrbreverse(nrbLeft);

        nrb1 = threePointCircle(point{24, iBarrier}, point{30, iBarrier}, point{29, iBarrier});
        nrb2 = threePointCircle(point{30, iBarrier}, point{31, iBarrier}, point{29, iBarrier});
        nrbRight = nrbglue(nrb1, nrb2, 2, 1);
        nrbRight.knots = nrbRight.knots / max(nrbRight.knots);
        % nrbRight = nrbreverse(nrbRight);

        srf(2, iBarrier) = nrbruled(nrbLeft, nrbRight);

        % flux barrier 3
        nrb1 = threePointCircle(point{27, iBarrier}, point{25, iBarrier}, point{24, iBarrier});
        nrb2 = threePointCircle(point{25, iBarrier}, point{5, iBarrier}, [0, 0]);
        nrb3 = threePointCircle(point{5, iBarrier}, point{7, iBarrier}, point{4, iBarrier});
        nrbRight = nrbglue(nrbglue(nrb1, nrb2, 2, 1), nrb3, 2, 1);
        nrbRight.knots = nrbRight.knots / max(nrbRight.knots);
        nrbLeft = nrbline(point{24, iBarrier}, point{4, iBarrier});
        srf(3, iBarrier) = nrbruled(nrbLeft, nrbRight);

        % flux barrier 4
        nrb1 = nrbextract(srf(1, iBarrier), 4);
        nrb2 = nrbextract(srf(2, iBarrier), 4);
        nrb3 = nrbextract(srf(3, iBarrier), 3);
        nrb4 = nrbline(point{11, iBarrier}, point{31, iBarrier});
        srf(4, iBarrier) = nrbcoons(nrb4, nrbreverse(nrb3), nrbreverse(nrb1), nrbreverse(nrb2));

        % flux barrier 5
        nrb1 = threePointCircle(point{14, iBarrier}, point{1, iBarrier}, point{3, iBarrier});
        nrb2 = threePointCircle(point{11, iBarrier}, point{8, iBarrier}, point{9, iBarrier});
        srf(5, iBarrier) = nrbruled(nrb1, nrb2);

        % flux barrier 6
        nrb1 = threePointCircle(point{11, iBarrier}, point{8, iBarrier}, point{9, iBarrier});
        nrb2 = threePointCircle(point{31, iBarrier}, point{28, iBarrier}, point{29, iBarrier});
        nrb3 = nrbextract(srf(4, iBarrier), 3);
        nrb4 = nrbline(point{8, iBarrier}, point{28, iBarrier});
        srf(6, iBarrier) = nrbcoons(nrb4, nrb3, nrbreverse(nrb1), nrbreverse(nrb2));

        % flux barrier 7
        nrb1 = threePointCircle(point{31, iBarrier}, point{28, iBarrier}, point{29, iBarrier});
        nrb2 = threePointCircle(point{34, iBarrier}, point{21, iBarrier}, point{23, iBarrier});
        srf(7, iBarrier) = nrbruled(nrb1, nrb2);

        % iron bridge 1
        nrb1 = threePointCircle(point{32, iBarrier}, point{12, iBarrier}, [0,0]);
        nrb2 = nrbextract(srf(3, iBarrier), 4);
        srf(8, iBarrier) = nrbruled(nrb1, nrb2);

        % air  gap 1
        nrb1 = threePointCircle(point{33, iBarrier}, point{13, iBarrier}, [0,0]);
        nrb2 = nrbextract(srf(8, iBarrier), 3);
        srf(9, iBarrier) = nrbruled(nrb1, nrb2);

        nBarrierPatches = 13;
        regions.Air = [regions.Air, nBarrierPatches*(iBarrier-1)+[1:7,9]];
        regions.Iron = [regions.Iron, nBarrierPatches*(iBarrier-1)+8];
        boundaries.Coupling.patches = [boundaries.Coupling.patches, nBarrierPatches*(iBarrier-1)+9];
        boundaries.Coupling.sides = [boundaries.Coupling.sides, 3];
    end

    for iBarrier = 1:params.BN
        if iBarrier < params.BN
            % flux carrier 1
            nrb1 = nrbextract(srf(7, iBarrier), 4); 
            nrb2 = nrbextract(srf(5, iBarrier+1), 3); 
            srf(10, iBarrier) = nrbruled(nrb1, nrb2);
            srf(10, iBarrier) = nrbdegelev(srf(10, iBarrier), [0, 1]);

            % flux carrier 2
            nrb1 = nrbextract(srf(2, iBarrier), 3); 
            nrb2 = nrbextract(srf(1, iBarrier+1), 3); 
            srf(11, iBarrier) = nrbruled(nrb1, nrb2);
            srf(11, iBarrier) = nrbdegelev(srf(11, iBarrier), [0, 1]);

            % iron bridge 2
            nrb1 = threePointCircle(point{12, iBarrier+1}, point{32, iBarrier}, [0,0]);
            nrb2 = nrbline(point{7, iBarrier+1}, point{27, iBarrier});
            srf(12, iBarrier) = nrbruled(nrb1, nrb2);

            % air gap 2
            nrb1 = threePointCircle(point{13, iBarrier+1}, point{33, iBarrier}, [0,0]);
            nrb2 = threePointCircle(point{12, iBarrier+1}, point{32, iBarrier}, [0,0]);
            srf(13, iBarrier) = nrbruled(nrb1, nrb2);

            regions.Air = [regions.Air, nBarrierPatches*(iBarrier-1)+13];
            regions.Iron = [regions.Iron, nBarrierPatches*(iBarrier-1)+[10,11,12]];
            boundaries.Coupling.patches = [boundaries.Coupling.patches, nBarrierPatches*(iBarrier-1)+13];
            boundaries.Coupling.sides = [boundaries.Coupling.sides, 3];
        end
    end

    % last center flux carrier
    a1 = atan2(point{26, iBarrier}(2), point{26, iBarrier}(1));
    % p1In = [cos(a1/2), sin(a1/2)]* (params.RRO - params.WOB(iBarrier));
    p1Mid = [cos(a1/2), sin(a1/2)]* params.RRO;
    p1Out = [cos(a1/2), sin(a1/2)]* params.RAG;
    % p2In = [params.RRO - params.WOB(iBarrier), 0];
    p2Mid = [params.RRO, 0];
    p2Out = [params.RAG, 0];

    nrb1 = threePointCircle(p1Mid, point{32, iBarrier}, [0,0]);
    nrb2 = nrbreverse(nrbextract(srf(2, iBarrier), 3));
    srf(10, iBarrier) = nrbruled(nrb1, nrb2);

    nrb1 = threePointCircle(p2Mid, p1Mid, [0,0]);
    nrb2 = nrbreverse(nrbextract(srf(7, iBarrier), 4));
    srf(11, iBarrier) = nrbruled(nrb1, nrb2);

    nrb1 = threePointCircle(p1Out, point{33, iBarrier}, [0,0]);
    nrb2 = nrbextract(srf(10, iBarrier), 3);
    srf(12, iBarrier) = nrbruled(nrb1, nrb2);

    nrb1 = threePointCircle(p2Out, p1Out, [0,0]);
    nrb2 = nrbextract(srf(11, iBarrier), 3);
    srf(13, iBarrier) = nrbruled(nrb1, nrb2);

    srf = reshape(srf, [] ,1);

    % inner iron part
    nrb1 = nrbcirc(params.RRI, [0,0], 0, beta);
    nrb2 = nrbreverse(nrbextract(srf(5), 3));
    nrb = nrbdegelev(nrbruled(nrb1, nrb2), [0, 1]);
    srf(end+1) = nrb;

    p1 = [cos(beta), sin(beta)]*params.RRI;
    p2 = [cos(beta), sin(beta)]*vecnorm(point{7,1});
    pm = 0.9*p2 + 0.1*p1;
    nrb1 = nrbglue(nrbdegelev(nrbline(p2, pm), 1), nrbdegelev(nrbline(pm, p1), 1), 2, 1);
    nrb1.knots = nrb1.knots/max(nrb1.knots);
    nrb2 = nrbextract(srf(1), 3);
    nrb = nrbdegelev(nrbruled(nrb1, nrb2), [0, 1]);
    srf(end+1) = nrb;

    p1 = [cos(beta), sin(beta)]*params.RRO;
    nrb1 = threePointCircle(point{12,1}, p1, [0,0]);
    nrb2 = nrbreverse(nrbextract(srf(end), 1));
    srf(end+1) = nrbruled(nrb1, nrb2);

    % last air part
    p1 = [cos(beta), sin(beta)]*params.RAG;
    nrb1 = threePointCircle(point{13,1}, p1, [0,0]);
    nrb2 = nrbextract(srf(end), 3);
    srf(end+1) = nrbruled(nrb1, nrb2);

    regions.Air = [regions.Air, nBarrierPatches*(iBarrier-1)+[12,13,17]];
    regions.Iron = [regions.Iron, nBarrierPatches*(iBarrier-1)+[10,11,14,15,16]];

    boundaries.Coupling.patches = [boundaries.Coupling.patches, nBarrierPatches*(iBarrier-1)+[12,13,17]];
    boundaries.Coupling.sides = [boundaries.Coupling.sides, [3,3,3]];
    boundaries.Dirichlet.patches = nBarrierPatches*(iBarrier-1)+14;
    boundaries.Dirichlet.sides = 3;

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
        nsrfHalf = nsrf;    % for second part of periodic BC
        nsrf = numel(srf);
        boundaries.Periodic1.patches = nsrf-[0,1,2];
        boundaries.Periodic1.sides = [2,2,3];
    end

    % rotate by half segment
    nsrf = numel(srf);
    for iSrf = 1:nsrf
        srf(iSrf) = nrbtform(srf(iSrf), vecrotz(beta));
    end

    % copy and rotate
    nsrf = numel(srf);
    regionsInit = regions;
    boundariesInit = boundaries;
    for I = 1:params.POLES/params.Symmetry-1
        for iSrf = 1:nsrf
            srf(I*nsrf+iSrf) = nrbtform(srf(iSrf), vecrotz(I*2*pi/params.POLES));
        end
        regions.Air = [regions.Air, regionsInit.Air+nsrf*I];
        regions.Iron = [regions.Iron, regionsInit.Iron+nsrf*I];
        boundaries.Dirichlet.patches = [boundaries.Dirichlet.patches, boundariesInit.Dirichlet.patches+nsrf*I];
        boundaries.Dirichlet.sides = [boundaries.Dirichlet.sides, boundariesInit.Dirichlet.sides];
        boundaries.Coupling.patches = [boundaries.Coupling.patches, boundariesInit.Coupling.patches+nsrf*I];
        boundaries.Coupling.sides = [boundaries.Coupling.sides, boundariesInit.Coupling.sides];
    end

    nsrf = numel(srf);
    if params.Symmetry ~= 1
        boundaries.Periodic2.patches = nsrf-nsrfHalf-[0, 1, 2];
        boundaries.Periodic2.sides = [2, 2, 3];
    end

    if (params.preview)
        preview_nurbs(srf, regions);
        % hold on
        % for iPoint = 1:numel(point)
        %     if isempty(point{iPoint})
        %         continue
        %     end
        %     scatter(point{iPoint}(1), point{iPoint}(2), "black", "filled");
        %     text(point{iPoint}(1)+1e-4, point{iPoint}(2), string(iPoint))
        % end
        % view(2)
        % axis equal
    end
end

% p1: startpoint, p2:endpoint, p3: center
function circ = threePointCircle(p1, p2, p3)
    angle1 = atan2(p1(2)-p3(2), p1(1)-p3(1));
    angle2 = atan2(p2(2)-p3(2), p2(1)-p3(1));
    R = vecnorm(p3-p1);

    assert(abs(vecnorm(p1-p3) - vecnorm(p2-p3)) < 1e-10);
    circ = nrbcirc(R, p3, angle1, angle2);
end