% https://www.mathcha.io/editor/dDK5WtP8twWsomy3WVUqXlB90CKBN6EqtKxQ9dP

function [srf, regions, boundaries, params] = GEO_DISstator_2(params)
    arguments
        params.preview (1,1) logical = true;
        params.Symmetry = 4;
        params.SLOTS = 24;              % Nuber of stator slots
        params.RAG = 50.5e-3;           % Radius Air Gap
        params.SRI = 51e-3;             % Stator Radius Inner
        params.SRO = 102e-3;            % Stator Radius Outer
        
        params.HSO = 1.6e-3;            % Height Slot Opening
        params.HY = 25.5e-3;            % Height Yoke
        params.WSO = 3.6e-3;            % Width Slot Opening
        params.WT = 6.0e-3;             % Width Tooth
        params.LumpWindings (1,1) logical = false;
    end
    
    % patch definitions
    regions.Iron = [2,3,4,6,8,12,24,20,18,16,15,14];
    regions.Air = [1,5,9,10,21,22,17,13];
    regions.Slot1up = [23,19];
    regions.Slot1down = [7,11];
    regions.Windings = [regions.Slot1up, regions.Slot1down];

    % calculated parameters
    a1 = asin((params.WT/2)/(params.SRI));
    a2 = asin((params.WT/2)/(params.SRI+params.HSO));
    a3 = asin((params.WT/2)/(params.SRO-params.HY));
    a4 = 2*pi/params.SLOTS - 2*a1;
    a5 = 2*pi/params.SLOTS - 2*a2;
    a6 = 2*pi/params.SLOTS - 2*a3;
    a7 = 2*asin(params.WSO/(2*params.SRI));
    a8 = 2*asin((params.WSO/2)/(params.SRI+params.HSO));

    srf(1) = nrbruled(nrbcirc(params.RAG,[0,0],0,a1),nrbcirc(params.SRI,[0,0],0,a1));
    srf(2) = nrbruled(nrbcirc(params.SRI,[0,0],0,a1),nrbcirc(params.SRI+params.HSO,[0,0],0,a2));
    srf(3) = nrbruled(nrbcirc(params.SRI+params.HSO,[0,0],0,a2),nrbcirc(params.SRO-params.HY,[0,0],0,a3));
    srf(4) = nrbruled(nrbcirc(params.SRO-params.HY,[0,0],0,a3),nrbcirc(params.SRO,[0,0],0,a3));

    srf(5) = nrbruled(nrbcirc(params.RAG,[0,0],a1,a1+(a4-a7)/2),nrbcirc(params.SRI,[0,0],a1,a1+(a4-a7)/2));
    srf(6) = nrbruled(nrbcirc(params.SRI,[0,0],a1,a1+(a4-a7)/2),nrbcirc(params.SRI+params.HSO,[0,0],a2,a2+(a5-a8)/2));
    srf(7) = nrbruled(nrbcirc(params.SRI+params.HSO,[0,0],a2,a2+(a5-a8)/2),nrbcirc(params.SRO-params.HY,[0,0],a3,a3+(a6-a8)/2));
    srf(8) = nrbruled(nrbcirc(params.SRO-params.HY,[0,0],a3,a3+(a6-a8)/2),nrbcirc(params.SRO,[0,0],a3,a3+(a6-a8)/2));

    srf(9) = nrbruled(nrbcirc(params.RAG,[0,0],a1+(a4-a7)/2,a1+(a4-a7)/2+a7/2),nrbcirc(params.SRI,[0,0],a1+(a4-a7)/2,a1+(a4-a7)/2+a7/2));
    srf(10) = nrbruled(nrbcirc(params.SRI,[0,0],a1+(a4-a7)/2,a1+(a4-a7)/2+a7/2),nrbcirc(params.SRI+params.HSO,[0,0],a2+(a5-a8)/2,a2+(a5-a8)/2+a8/2));
    srf(11) = nrbruled(nrbcirc(params.SRI+params.HSO,[0,0],a2+(a5-a8)/2,a2+(a5-a8)/2+a8/2),nrbcirc(params.SRO-params.HY,[0,0],a3+(a6-a8)/2,a3+(a6-a8)/2+a8/2));
    srf(12) = nrbruled(nrbcirc(params.SRO-params.HY,[0,0],a3+(a6-a8)/2,a3+(a6-a8)/2+a8/2),nrbcirc(params.SRO,[0,0],a3+(a6-a8)/2,a3+(a6-a8)/2+a8/2));
    % mirror first half winding
    nsrf = numel(srf);
    mirror = diag(ones(4,1));
    mirror(2,2) = -1;
    for i = 1:nsrf
        srf(nsrf+i) = nrbtform(srf(i), mirror*vecrotz(-(2*pi/params.SLOTS)));
    end

    % Boundaries
    boundaries.Dirichlet.patches = [4,8,12,24,20,16];
    boundaries.Dirichlet.sides = [4,4,4,4,4,4];
    boundaries.Coupling.patches = [1,5,9,21,17,13];
    boundaries.Coupling.sides = [3,3,3,3,3,3];
    
    % Copy and rotate
    nsrf = numel(srf);
    IronInit = regions.Iron;
    AirInit = regions.Air;
    WindingsInit = regions.Windings;
    boundariesInit = boundaries;
    for I = 1:params.SLOTS/params.Symmetry-1
        for i = 1:nsrf
            srf(I*nsrf+i) = nrbtform(srf(i), vecrotz(I*2*pi/params.SLOTS));
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
    if params.Symmetry == 1
        % do nothing
    else
        boundaries.Periodic1.patches = [1,2,3,4];
        boundaries.Periodic1.sides = [1,1,1,1];
        boundaries.Periodic2.patches = nsrf-[11,10,9,8];
        boundaries.Periodic2.sides = [1,1,1,1];
    end

    if (params.preview)
        preview_nurbs(srf, regions);
    end
end