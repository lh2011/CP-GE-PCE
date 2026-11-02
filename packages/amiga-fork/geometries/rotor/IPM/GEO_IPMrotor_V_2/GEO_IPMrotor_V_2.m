

function [srf, regions, boundaries, params, magnetization] = GEO_IPMrotor_V_2 (params)
    arguments
        params.preview = true;
        params.POLES = 8;
        params.Symmetry = 8;
        params.REFINEMENT = [0, 0, 3, 1];
    
        params.RRO = 78.6e-3;           % Rotor Radius Outer
        params.RRI = 26.5e-3;           % Rotor Radius Inner
        params.RAG = 78.9e-3;           % Radius Air Gap
    
        params.DMAG = 70.5e-3;          % Distance Magnet Positioning 
        params.LSLIT1 = 1.9e-3;         % Air Slit Positioning
        params.LSLIT2 = 6e-3;           % Air Slit Positioning
        params.DSLIT1 = 3.5e-3;         % Air Slit Positioning
        params.DSLIT2 = 0.1e-3;         % Air Slit Positioning
        params.MA = 160;                % Magnet opening Angle
        params.MH = 4.0e-3;             % Magnet Height
        params.MW = 21.0e-3;            % Magnet Width
        params.RA1 = 190;               % Air Slit Angle 1
        params.RA2 = 210.0;             % Air Slit Angle 2
        params.RS = 1.4e-3;             % Rotor Iron Slit thickness
        params.WMAG = 1.7e-3;           % Width Magnet Positoning
    end

    % calculated parameters
    motorAngle = 2*pi/params.POLES;

    % magnet
    point{1} = [params.DMAG, params.WMAG];
    point{2} = point{1} + params.MH*[-cosd(90-params.MA/2), sind(90-params.MA/2)];
    point{3} = point{1} + params.MW*[cosd(params.MA/2), sind(params.MA/2)];
    point{4} = point{3} + params.MH*[-cosd(90-params.MA/2), sind(90-params.MA/2)];
    % inner air slit
    point{5} = [params.DMAG-params.DSLIT2, params.RS/2];
    point{6} = [params.DMAG-params.DSLIT1, params.RS/2];
    point{7} = [point{5}(1), 0];
    point{8} = [point{6}(1), 0];
    % outer air slit
    point{9} = point{3} + params.LSLIT2*[cosd(params.RA2-180+params.MA/2), sind(params.RA2-180+params.MA/2)];
    point{10} = point{4} + params.LSLIT1*[cosd(params.RA2-180+params.MA/2), sind(params.RA2-180+params.MA/2)];
    point{11} = (params.RRI + params.RRO)/2*[cos(motorAngle/2), sin(motorAngle/2)];

    % hard coded (same discretization independent of parameters) / vs projected
    % a1 = motorAngle/8; 
    % a2 = 2*motorAngle/8;  
    % a3 = 3*motorAngle/8; 
    % a4 = motorAngle/2;
    a1 = atan2(point{5}(2), point{5}(1));
    a2 = atan2(point{1}(2), point{1}(1));
    a3 = atan2(point{3}(2), point{3}(1));
    a4 = motorAngle/2;

    % circs defining the refinement of the motor
    circ1 = nrbcirc(params.RRO, [0, 0], 0, a1);
    circ2 = nrbcirc(params.RRO, [0, 0], a1, a2);
    circ3 = nrbcirc(params.RRO, [0, 0], a2, a3);
    circ4 = nrbcirc(params.RRO, [0, 0], a3, a4);
    % refine circle as wanted...
    circ1 = nrbrefine1D(circ1, params.REFINEMENT(1));
    circ2 = nrbrefine1D(circ2, params.REFINEMENT(2));
    circ3 = nrbrefine1D(circ3, params.REFINEMENT(3));
    circ4 = nrbrefine1D(circ4, params.REFINEMENT(4));
% 
%     % Hard coded values for air gap, which needs to remain the same
    A1 = motorAngle/20; 
    A2 = motorAngle/10;
    A3 = motorAngle/3; 
    A4 = motorAngle/2;

    % Air pacthces
    srf(1) = nrbruled(nrbcirc(params.RAG, [0,0], 0,  A1), circ1);
    srf(2) = nrbruled(nrbcirc(params.RAG, [0,0], A1, A2), circ2);
    srf(3) = nrbruled(nrbcirc(params.RAG, [0,0], A2, A3), circ3);
    srf(4) = nrbruled(nrbcirc(params.RAG, [0,0], A3, A4), circ4);

    srf(5) = nrbruled(nrbextract(srf(1), 4), nrbline(point{7}, point{5}));
    srf(6) = nrbruled(nrbextract(srf(2), 4), nrbline(point{5}, point{1}));
    srf(7) = nrbruled(nrbextract(srf(3), 4), nrbline(point{1}, point{3}));
    srf(8) = nrbruled(nrbextract(srf(4), 4), nrbline(point{3}, point{9}));
    
    srf(9) = nrbruled(nrbextract(srf(5), 4), nrbline(point{8}, point{6}));
    srf(10) = nrbruled(nrbextract(srf(6), 4), nrbline(point{6}, point{2}));
    srf(11) = nrbruled(nrbextract(srf(7), 4), nrbline(point{2}, point{4}));
    srf(12) = nrbruled(nrbextract(srf(8), 4), nrbline(point{4}, point{10}));
    srf(13) = nrbruled(nrbreverse(nrbextract(srf(8), 2)), nrbline(point{10}, point{11}));

    srf(14) = nrbruled(nrbextract(srf(9), 4), nrbcirc(params.RRI, [0,0], 0, a1));
    srf(15) = nrbruled(nrbextract(srf(10), 4), nrbcirc(params.RRI, [0,0], a1, (a2+a3)/2));
    srf(16) = nrbruled(nrbextract(srf(11), 4), nrbcirc(params.RRI, [0,0], (a2+a3)/2, a4));
    srf(17) = nrbruled(nrbextract(srf(12), 4), nrbline(params.RRI*[cos(motorAngle/2), sin(motorAngle/2)], point{11}));

    % mirror patches
    nsrf = numel(srf);
    mirror = diag(ones(4, 1));
    mirror(2, 2) = -1;
    for iSrf = 1:nsrf
        srf(nsrf+iSrf) = nrbtform(srf(iSrf), mirror);
    end
    % Rotate again by half a slot and upward facing
    nsrf = numel(srf);
    for iSrf = 1:nsrf
        srf(iSrf) = nrbtform(srf(iSrf), vecrotz(pi/params.POLES));
    end

    % region definition
    regions.Iron = [5:9,13:17,22:26,30:34];
    regions.Magnet1down = 28;
    regions.Magnet1up = 11;
    % regions.Magnets = [regions.Magnet1down, regions.Magnet1up];
    regions.Air = [1:4,18:21, 10,12,27,29];

    % boundary definition
    boundaries.Dirichlet.patches = [14,15,16,31,32,33];
    boundaries.Dirichlet.sides = [4,4,4,4,4,4];
    boundaries.Coupling.patches = [1,2,3,4,18,19,20,21];
    boundaries.Coupling.sides = [3,3,3,3,3,3,3,3];
    % Magnetization
    magnetization = [pi/2-deg2rad(params.MA/2)+pi/params.POLES, -pi/2+deg2rad(params.MA/2)+pi/params.POLES];

    % Copy and rotate
    nsrf = numel(srf);
    IronInit = regions.Iron;
    AirInit = regions.Air;
    % MagnetsInit = regions.Magnets;
    boundariesInit = boundaries;
    magnetizationInit = magnetization;
    for I = 1:params.POLES/params.Symmetry-1
        for i = 1:nsrf
            srf(I*nsrf+i) = nrbtform(srf(i), vecrotz(I*2*pi/params.POLES));
        end
        regions.Iron = [regions.Iron, IronInit+nsrf*I];
        regions.Air = [regions.Air, AirInit+nsrf*I];
        regions.("Magnet"+num2str(I+1)+"down") = regions.Magnet1down + nsrf*I;
        regions.("Magnet"+num2str(I+1)+"up") = regions.Magnet1up + nsrf*I;
        % regions.Magnets = [regions.Magnets, MagnetsInit+ nsrf*I];
        boundaries.Dirichlet.patches = [boundaries.Dirichlet.patches, boundariesInit.Dirichlet.patches+nsrf*I];
        boundaries.Dirichlet.sides = [boundaries.Dirichlet.sides, boundariesInit.Dirichlet.sides];
        boundaries.Coupling.patches = [boundaries.Coupling.patches, boundariesInit.Coupling.patches+nsrf*I];
        boundaries.Coupling.sides = [boundaries.Coupling.sides, boundariesInit.Coupling.sides];
        magnetization = [magnetization, magnetizationInit+(pi+motorAngle)*I];
    end

    % Additional Boundaries if not full machine is simulated
    if params.Symmetry ~= 1
        nsrf = numel(srf);
        boundaries.Periodic1.patches = [21,30,34];
        boundaries.Periodic1.sides = [2,2,4];
        boundaries.Periodic2.patches = nsrf-[30,21,17];
        boundaries.Periodic2.sides = [2,2,4];
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