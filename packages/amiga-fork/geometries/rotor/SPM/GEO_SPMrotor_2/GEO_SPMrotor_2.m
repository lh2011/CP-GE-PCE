

function [srf, regions, boundaries, params] = GEO_SPMrotor_2(params)
    arguments
        params.preview = true;
        params.POLES = 10;
        params.Symmetry = 2;
    
        params.HM = 5.3e-3;             % Height of Magnet in Sum
        params.HM1 = 0.6e-3;            % Height of Magnet Part Outer
        params.HM2 = 1.2e-3;            % Height of Magnet Part Inner
        params.RHY = 3.5e-3;            % Rotor Height Yoke
        params.RRI = 43e-3;             % Rotor Radius Inner
        params.RRO = [];                % Rotor Radius Outer, will be calculated 
        params.RAG = 42.5e-3;           % Radius Air Gap
        params.RatioHM = 0.95;          % Magnet Ratio of circumference
        params.RatioHM1 = 0.9;          % Magnet Ratio of circumference
        params.RatioHM2 = 0.85;          % Magnet Ratio of circumference
    end
    params.RRO = params.RRI + params.HM + params.RHY;

    center = [0,0,0];
    beta = pi/params.POLES;
    % Inner layer Air
    srf(1) = nrbruled(nrbcirc(params.RAG, center, 0, beta*params.RatioHM2), nrbcirc(params.RRI, center, 0, beta*params.RatioHM2));
    srf(2) = nrbruled(nrbcirc(params.RAG, center, beta*params.RatioHM2, beta), nrbcirc(params.RRI, center, beta*params.RatioHM2, beta));
    % Second layer Magnet / Air
    srf(3) = nrbruled(nrbcirc(params.RRI, center, 0, beta*params.RatioHM2), nrbcirc(params.RRI+params.HM2, center, 0, beta*params.RatioHM));
    srf(4) = nrbruled(nrbcirc(params.RRI, center, beta*params.RatioHM2, beta), nrbcirc(params.RRI+params.HM2, center, beta*params.RatioHM, beta));
    % Third layer Magnet / Air
    srf(5) = nrbruled(nrbcirc(params.RRI+params.HM2, center, 0, beta*params.RatioHM), nrbcirc(params.RRI + params.HM - params.HM1, center, 0, beta*params.RatioHM));
    srf(6) = nrbruled(nrbcirc(params.RRI+params.HM2, center, beta*params.RatioHM, beta), nrbcirc(params.RRI + params.HM - params.HM1, center, beta*params.RatioHM, beta));
    % Forth layer Magnet / Air
    srf(7) = nrbruled(nrbcirc(params.RRI + params.HM - params.HM1, center, 0, beta*params.RatioHM), nrbcirc(params.RRI + params.HM, center, 0, beta*params.RatioHM1));
    srf(8) = nrbruled(nrbcirc(params.RRI + params.HM - params.HM1, center, beta*params.RatioHM, beta), nrbcirc(params.RRI + params.HM, center, beta*params.RatioHM1, beta));
    % Outer layer Iron
    srf(9) = nrbruled(nrbcirc(params.RRI + params.HM, center, 0, beta*params.RatioHM1), nrbcirc(params.RRI + params.HM + params.RHY, center, 0, beta*params.RatioHM1));
    srf(10) = nrbruled(nrbcirc(params.RRI + params.HM , center, beta*params.RatioHM1, beta), nrbcirc(params.RRI + params.HM + params.RHY, center, beta*params.RatioHM1, beta));

    % mirror patches
    nsrf = numel(srf);
    mirror = diag(ones(4, 1));
    mirror(2, 2) = -1;
    for iSrf = 1:nsrf
        srf(nsrf+iSrf) = nrbtform(srf(iSrf), mirror);
    end

    % rotate by half segment
    nsrf = numel(srf);
    for iSrf = 1:nsrf
        srf(iSrf) = nrbtform(srf(iSrf), vecrotz(beta));
    end

    % region definition
    regions.Iron = [9,10,19,20];
    regions.Air = [1,2,4,6,8, 11,12,14,16,18];
    regions.Magnet1 = [3,5,7,13,15,17];

    % boundary definition
    boundaries.Coupling.patches = [1,2,11,12];
    boundaries.Coupling.sides = [3,3,3,3];
    boundaries.Dirichlet.patches = [9,10,19,20];
    boundaries.Dirichlet.sides = [4,4,4,4];
    
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
        regions.("Magnet"+num2str(I+1)) = regions.Magnet1 + nsrf*I;
        boundaries.Dirichlet.patches = [boundaries.Dirichlet.patches, boundariesInit.Dirichlet.patches+nsrf*I];
        boundaries.Dirichlet.sides = [boundaries.Dirichlet.sides, boundariesInit.Dirichlet.sides];
        boundaries.Coupling.patches = [boundaries.Coupling.patches, boundariesInit.Coupling.patches+nsrf*I];
        boundaries.Coupling.sides = [boundaries.Coupling.sides, boundariesInit.Coupling.sides];
    end

    nsrf = numel(srf);
    if params.Symmetry ~= 1
        boundaries.Periodic1.patches = [12,14,16,18,20];
        boundaries.Periodic1.sides = [2,2,2,2,2];
        boundaries.Periodic2.patches = nsrf-[18,16,14,12,10];
        boundaries.Periodic2.sides = [2,2,2,2,2];
    end
    if (params.preview)
        preview_nurbs(srf, regions);
    end
end