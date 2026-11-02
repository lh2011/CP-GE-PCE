% https://www.mathcha.io/editor/QJWQwfM8fqNU2r1qweIG4MjovtpX8G3Ds24pk7Y

function [srf, regions, boundaries, params] = GEO_SPMrotor_1(params)
    arguments
        params.preview (1,1) logical = true;
        params.POLES = 4;
        params.Symmetry = 2;
    
        params.HM = 4.35e-3;            % Height of Magnet
        params.RRI = 8e-3;              % Rotor Radius Inner
        params.RRO = 47e-3/2;           % Rotor Radius Outer
        params.RAG = 23.7e-3;           % Radius Air Gap
        params.MR = 0.5;                % Magnet Ratio of circumference
    end

    center = [0,0,0];
    beta = pi/params.POLES;
    srf(1) = nrbruled(nrbcirc(params.RRI, center, 0, beta*params.MR), nrbcirc(params.RRO-params.HM, center, 0, beta*params.MR));
    srf(2) = nrbruled(nrbcirc(params.RRO-params.HM, center, 0, beta*params.MR), nrbcirc(params.RRO, center, 0, beta*params.MR));
    srf(3) = nrbruled(nrbcirc(params.RRO, center, 0, beta*params.MR), nrbcirc(params.RAG, center, 0, beta*params.MR));
    
    srf(4) = nrbruled(nrbcirc(params.RRI, center, beta*params.MR, beta), nrbcirc(params.RRO-params.HM, center, beta*params.MR, beta));
    srf(5) = nrbruled(nrbcirc(params.RRO-params.HM, center, beta*params.MR, beta), nrbcirc(params.RRO, center, beta*params.MR, beta));
    srf(6) = nrbruled(nrbcirc(params.RRO, center, beta*params.MR, beta), nrbcirc(params.RAG, center, beta*params.MR, beta));
    
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
    regions.Iron = [1,4,5,7,10,11];
    regions.Air = [3,6,9,12];
    regions.Magnet1 = [2,8];

    % boundary definition
    boundaries.Dirichlet.patches = [1,4,7,10];
    boundaries.Dirichlet.sides = [3,3,3,3];
    boundaries.Coupling.patches = [3,6,9,12];
    boundaries.Coupling.sides = [4,4,4,4];
    
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
        boundaries.Periodic1.patches = [10,11,12];
        boundaries.Periodic1.sides = [2,2,2];
        boundaries.Periodic2.patches = nsrf-[8,7,6];
        boundaries.Periodic2.sides = [2,2,2];
    end
    
    if (params.preview)
        preview_nurbs(srf, regions);
    end
end