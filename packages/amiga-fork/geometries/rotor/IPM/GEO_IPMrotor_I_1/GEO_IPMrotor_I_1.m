% https://www.mathcha.io/editor/P1GLZsWPU5OtxQ7E6rh1q5B6PTDw7LEpTExPn4G

function [srf, regions, boundaries, params, magnetization] = GEO_IPMrotor_I_1(params)
    arguments
        params.preview (1,1) logical = true;
        params.POLES       = 6;                 % Number of poles
        params.Symmetry    = 6;                 % Reduction due to symmetry
        params.MW          = 19e-3;             % Magnet Width
        params.MH          = 7e-3;              % Magnet Height
        params.MAG         = 7e-3;              % Magnet to Air Gap distance
        params.RRO         = 44e-3;             % Rotor Radius Outer
        params.RRI         = 16e-3;             % Rotor Radius Inner
        params.RAG         = 44.5e-3;           % Air Gap Radius
        params.Theta1      = 0.37;              % slit positioning
        params.Theta2      = 0.00558;           % slit positioning
        params.LW1         = 2e-3;              % slit positioning
        params.LW2         = 3e-3;              % slit positioning
        params.REFINEMENT  = [1,2,3];           % increases number of control points
        params.Type {mustBeMember(params.Type, ["linear", "nonlinear"])}       = "nonlinear";   % material definition
    end

    % Derived parameters
    beta = 2*pi/params.POLES;
    alpha = beta/2;
    center = [0, 0, 0];
    dMag = params.RRO - params.MH - params.MAG;
    r6 = params.RRO - params.LW2;
    r7 = params.RRO - params.LW1;
    % fixing the air gap for now
    THETA1 = 0.35*6/params.POLES;
    THETA2 = 0.4*6/params.POLES;
    Theta2 = params.Theta1+params.Theta2;
    
    % setting points
    point{1} = params.RRI*[cos(alpha), -sin(alpha)];
    point{2} = params.RRI*[cos(alpha), sin(alpha)];
    
    point{3} = [dMag, -params.MW/2];
    point{4} = [dMag, params.MW/2];
    
    point{5} = r6*[cos(alpha), -sin(alpha)];
    point{6} = r6*[cos(Theta2), -sin(Theta2)];
    point{7} = r7*[cos(params.Theta1), -sin(params.Theta1)];

    point{8} = [dMag + params.MH, -params.MW/2];
    point{9} = [dMag + params.MH, params.MW/2];

    point{10} = r7*[cos(params.Theta1), sin(params.Theta1)];    
    point{11} = r6*[cos(Theta2), -sin(Theta2)];
    point{12} = r6*[cos(alpha), -sin(alpha)];

    point{13} = params.RRO*[cos(alpha), -sin(alpha)];
    point{14} = params.RRO*[cos(Theta2), -sin(Theta2)];
    point{15} = params.RRO*[cos(params.Theta1), -sin(params.Theta1)];
    point{16} = params.RRO*[cos(params.Theta1), sin(params.Theta1)];
    point{17} = params.RRO*[cos(Theta2), sin(Theta2)];
    point{18} = params.RRO*[cos(alpha), sin(alpha)];
    
    point{19} = params.RAG*[cos(alpha), -sin(alpha)];
    point{20} = params.RAG*[cos(THETA2), -sin(THETA2)];
    point{21} = params.RAG*[cos(THETA1), -sin(THETA1)];
    point{22} = params.RAG*[cos(THETA1), sin(THETA1)];
    point{23} = params.RAG*[cos(THETA2), sin(THETA2)];
    point{24} = params.RAG*[cos(alpha), sin(alpha)];
    
    % Patch 1 - Bottom Right Iron Part
    LeftSide = nrbline(point{5}, point{6});
    RightSide  = nrbrefine1D(nrbcirc (params.RRO, center, -alpha, -Theta2), params.REFINEMENT(1));
    srf(1) = nrbruled(LeftSide, RightSide);

    % Patch 2 - Bottom Middle Air Slit
    LeftSide = nrbline(point{6}, point{7});
    RightSide = nrbcirc(params.RRO, center, -Theta2, -params.Theta1);
    srf(2) = nrbruled(LeftSide, RightSide);
    
    % Patch 3 - Right Iron Part
    LeftSide = nrbline(point{7}, point{8});
    RightSide  = nrbrefine1D(nrbcirc(params.RRO, center, -params.Theta1, -params.Theta1/2), params.REFINEMENT(2));
    srf(3) = nrbruled(LeftSide, RightSide);

    % Patch 4 - Center Iron Part
    LeftSide = nrbline(point{8}, point{9});
    RightSide  = nrbrefine1D(nrbcirc(params.RRO, center, -params.Theta1/2, params.Theta1/2), params.REFINEMENT(3));
    srf(4) = nrbruled(LeftSide, RightSide);

    % Patch 5 - Upper Iron Part
    mirrorx = diag([1,-1,1,1]);
    srf(5) = nrbtform(srf(3), mirrorx);
    
    % Patch 6 - Top Middle Air Slit
    srf(6) = nrbtform(srf(2), mirrorx);

    % Patch 7 - Top Right Iron Part
    srf(7) = nrbtform(srf(1), mirrorx);

    % Patch 8 - Bottom Left Air Slit
    LeftSide = nrbline(point{6}, point{3});
    TopSide = nrbline(point{3}, point{8});
    RightSide = nrbextract(srf(3), 3); %nrbline(point{7}, point{8});
    BotSide = nrbextract(srf(2), 3);  %nrbline(point{6}, point{7});
    srf(8) = nrbcoons(BotSide, TopSide, LeftSide, RightSide);

    % Patch 9 - Permanent Magnet
    LeftSide = nrbline(point{3}, point{4});
    TopSide = nrbline(point{4}, point{9});
    RightSide = nrbextract(srf(4), 3); %nrbline(point{8}, point{9});
    BotSide = nrbextract(srf(8), 4); %nrbline(point{3}, point{8});
    srf(9) = nrbcoons(BotSide, TopSide, LeftSide, RightSide);

    % Patch 10 - Bottom Left Iron Part
    LeftSide = nrbline(point{1}, point{3});
    TopSide = nrbreverse(nrbextract(srf(8), 1));
    RightSide = nrbextract(srf(1), 3); %nrbline(point{5}, point{6});
    BotSide = nrbline(point{1}, point{5});
    srf(10) = nrbcoons(BotSide, TopSide, LeftSide, RightSide);
    
    % Patch 11 - Left Iron Part
    TopSide = nrbline(point{2}, point{4});
    LeftSide  = nrbcirc(params.RRI, center, -alpha, alpha);
    RightSide = nrbextract(srf(9), 1); %nrbline(point{3}, point{4});
    BotSide = nrbextract(srf(10), 1);
    srf(11) = nrbcoons(BotSide, TopSide, LeftSide, RightSide);

    % Patch 12 - Top Left Iron Part
    srf(12) = nrbtform (srf(10), mirrorx);
    
    % Patch 13 - Top Left Air Slit
    srf(13) = nrbtform(srf(8), mirrorx);

    % Patch 14 - Bottom Right Air Slit
    LeftSide = nrbextract(srf(1), 4); 
    RightSide  = nrbcirc(params.RAG, center, -alpha, -THETA2);
    srf(14) = nrbruled(LeftSide, RightSide);

    % Patch 15 - Bottom Right Air Slit 2
    LeftSide = nrbextract(srf(2), 4); 
    RightSide  = nrbcirc(params.RAG, center, -THETA2, -THETA1);
    srf(15) = nrbruled(LeftSide, RightSide);

    % Patch 16 - Middle Rigth Air Slit 1
    LeftSide = nrbextract(srf(3), 4);
    RightSide  = nrbcirc(params.RAG, center, -THETA1, -THETA1/2);
    srf(16) = nrbruled(LeftSide, RightSide);

    % Patch 17 - Middle Rigth Air Slit 2
    LeftSide = nrbextract(srf(4), 4);
    RightSide  = nrbcirc(params.RAG, center, -THETA1/2, THETA1/2);
    srf(17) = nrbruled(LeftSide, RightSide);

    % Patch 18 - Middle Rigth Air Slit 3
    srf(18) = nrbtform(srf(16), mirrorx);

    % Patch 18 - Top Right Air Slit 1
    srf(19) = nrbtform(srf(15), mirrorx);

    % Patch 18 - Top Right Air Slit 2
    srf(20) = nrbtform(srf(14), mirrorx);

    nsrf = numel(srf);
    for iSrf = 1:nsrf
        srf(iSrf) = nrbtform(srf(iSrf), vecrotz(pi/params.POLES));
    end

    % region definition
    if params.Type == "linear"
        regions.Iron = [1,3,4,5,7,10,11,12];
        regions.Air = [2,6,8,13,14:20];
    elseif params.Type == "nonlinear"
        regions.Iron = [1,2,3,4,5,6,7,10,11,12];
        regions.Air = [8,13,14:20];
    end
    regions.Magnet1 = [9];

    % boundary definition
    boundaries.Dirichlet.patches = 11;
    boundaries.Dirichlet.sides = 1;
    boundaries.Coupling.patches = [14,15,16,17,18,19,20];
    boundaries.Coupling.sides = [4,4,4,4,4,4,4];
    % magnetization
    magnetization = beta/2;

    % Copy and rotate
    nsrf = numel(srf);
    regionsInit = regions;
    boundariesInit = boundaries;
    magnetizationInit = magnetization;
    for I = 1:params.POLES/params.Symmetry-1
        for iSrf = 1:nsrf
            srf(I*nsrf+iSrf) = nrbtform(srf(iSrf), vecrotz(I*2*pi/params.POLES));
        end
        regions.Iron = [regions.Iron, regionsInit.Iron+nsrf*I];
        regions.Air = [regions.Air, regionsInit.Air+nsrf*I];
        regions.("Magnet"+num2str(I+1)) = regions.Magnet1 + nsrf*I;
        boundaries.Dirichlet.patches = [boundaries.Dirichlet.patches, boundariesInit.Dirichlet.patches+nsrf*I];
        boundaries.Dirichlet.sides = [boundaries.Dirichlet.sides, boundariesInit.Dirichlet.sides];
        boundaries.Coupling.patches = [boundaries.Coupling.patches, boundariesInit.Coupling.patches+nsrf*I];
        boundaries.Coupling.sides = [boundaries.Coupling.sides, boundariesInit.Coupling.sides];
        magnetization = [magnetization, magnetizationInit+(pi+beta)*I];
    end
    % Rotate again by half a slot and upward facing
    % nsrf = numel(srf);
    % for I = 1:nsrf
    %     srf(I) = nrbtform(srf(I), vecrotz(pi/params.POLES + pi/2 - pi/params.Symmetry));
    % end

    nsrf = numel(srf);
    if params.Symmetry ~= 1
        boundaries.Periodic1.patches = [10,1,14];
        boundaries.Periodic1.sides = [3,1,1];
        boundaries.Periodic2.patches = nsrf-[8,13,0];
        boundaries.Periodic2.sides = [3,1,1];
    end

    if (params.preview)
        preview_nurbs(srf, regions);
        % for iSrf = 1:numel(point)
        %     scatter(point{iSrf}(1), point{iSrf}(2), "black", "filled");
        %     text(point{iSrf}(1)+1e-5, point{iSrf}(2), string(iSrf))
        % end
    end
end
