

function [srf, regions, boundaries, params] = GEO_HPWstator_1(params)
    arguments
        params.preview = true;
        params.Symmetry = 4;
        params.SLOTS = 36;             % Nuber of stator slots
        
        params.HSO = 1e-3;             % Height Slot Opening
        params.RAG = 50.5e-3;          % Radius Air Gap
        params.SRI = 51e-3;            % Stator Radius Inner
        params.SRO = 90e-3;            % Stator Radius Outer
        params.SR1 = 0.6e-3;           % Slot Radius
        params.WSO = 2.5e-3;           % Width Slot Opening
        
        params.PH = 4e-3;              % Pin Height
        params.PN = 5;                 % Number of Hairpins
        params.PP1 = 2e-4;             % Pin Padding between Pins
        params.PP2 = 2e-4;             % Pin Padding between Iron and Air
        params.PW = 5e-3;              % Pin Width
        params.PR = 1e-3;              % Pin Radius
        params.LumpWindings (1,1) logical = false;
    end
    
    % Overwrite standard parameters with input parameters
    if params.PN<2
        error("At least 2 Pins are required")
    end
    
    % calculated parameters
    center = [0, 0, 0];    
    beta = 2*pi/params.SLOTS;
    alpha = beta/2;
    gamma1 = asin(params.WSO/2 / params.SRI);
    Rslot = sqrt((params.SRI*cos(gamma1) + params.HSO)^2 + params.WSO^2/4);
    gamma2 = asin(params.WSO/2 / Rslot); 

    % region definitions
    regions.Iron = [4, 6, 21, 22];
    regions.Air = [1:3, 5, 7:9, 14,15,18];
    regions.Windings = [10:13, 16,17,19, 20];
    regions.Slot1_Pin1 = regions.Windings;
    
    %setting points
    point{1} = [params.SRI*cos(alpha), -params.SRI*sin(alpha)];
    point{2} = [params.SRI*cos(gamma1), -params.SRI*sin(gamma1)];
    point{3} = [params.SRI*cos(gamma1), params.SRI*sin(gamma1)];
    point{4} = [params.SRI*cos(alpha), params.SRI*sin(alpha)];

    point{5} = point{2} + [params.HSO, 0];
    point{6} = point{3} + [params.HSO, 0];
    point{7} = [point{5}(1) + params.SR1, -params.PW/2 - params.PP2 + params.SR1];
    point{8} = point{7} - params.SR1/sqrt(2);
    point{9} = point{7} - [params.SR1, 0];
    point{10} = [point{5}(1) +  params.PP1 + params.PR, -params.PW/2 + params.PR];
    point{11} = point{10} - params.PR/sqrt(2);
    point{12} = point{10} - [params.PR, 0];
    point{13} = point{5} + [params.PP1, 0];
    point{14} = [point{13}(1), -point{13}(2)]; % mirror
    
    fixedRatio = 1/3;   % ratio of outer lengths w.r.t inner rectangle length
    wa = 0.5*params.PW*fixedRatio;
    wb = params.PW - 2*wa;
    ha = 0.5*params.PH*fixedRatio;
    hb = params.PH - 2*ha;
    point{15} = [point{13}(1) + ha, -wb/2];
    point{16} = [point{15}(1), -point{15}(2)]; % mirror
    point{17} = [(point{13}(1) + point{15}(1))/2, point{13}(2)];
    point{18} = [point{17}(1), -point{17}(2)]; % mirror
    point{19} = point{7} - [0, params.SR1];
    point{20} = [point{5}(1) + params.PH + 1.5*params.PP1, -params.PW/2 - params.PP2];
    point{21} = point{10} - [0, params.PR];
    point{22} = point{21} + [params.PH - 2*params.PR, 0];
    point{23} = point{22} + [0, params.PR];
    point{24} = point{23} + [params.PR, 0];
    point{25} = [point{24}(1) - ha, point{15}(2)];
    point{26} = [point{25}(1), -point{25}(2)]; % mirror
    point{27} = [point{24}(1), -point{24}(2)]; % mirror
    point{28} = [point{23}(1), -point{23}(2)]; % mirror
    point{29} = [point{20}(1), -point{20}(2)]; % mirror
    point{30} = [point{24}(1), point{20}(2)];
    point{31} = [point{20}(1), point{24}(2)];
    point{32} = [point{31}(1), -point{31}(2)]; % mirror
    point{33} = point{29} - [params.PH + params.PP1, 0];
    point{34} = point{33} + [params.PP1/2, 0];
    point{35} = [point{30}(1), point{29}(2)];
    point{36} = point{20} + [params.PP1/2, 0];
    point{37} = point{29} + [params.PP1/2, 0];
    point{38} = point{36} + [0, params.SR1];
    point{39} = point{37} + [0, -params.SR1];
    
    point{40} = point{36} + [-params.SR1, params.SR1];
    point{41} = point{37} + [-params.SR1, -params.SR1];

    for iPin = 1:params.PN-1
        l = point{5}(1) + 0.5*params.PP1 + iPin* (params.PP1 + params.PH);
        point{41+iPin} = [l, l*tan(alpha)];
    end
    point{end+1} = params.SRO * [cos(alpha), sin(alpha)];

    
    % creating surfaces
    % Air inside COUPLING INTERFACE IS HARD CODED FOR SAME BOUNDARY REPRESENTATION
    srf(1) = nrbruled(nrbcirc(params.RAG, center, -beta/2, -beta/6), nrbcirc(params.SRI, center, -beta/2, -gamma1));
    srf(1) = nrbkntins(srf(1), {[0.5,0.5], []}); % adjust for same mapping as following patches
    srf(2) = nrbruled(nrbcirc(params.RAG, center, -beta/6, beta/6), nrbcirc(params.SRI, center, -gamma1, gamma1));
    srf(2) = nrbkntins(srf(2), {[1/3, 1/3, 2/3, 2/3], []});
    srf(3) = nrbruled(nrbcirc(params.RAG, center, beta/6, beta/2), nrbcirc(params.SRI, center, gamma1, beta/2));
    srf(3) = nrbkntins(srf(3), {[0.5,0.5], []}); % adjust for same mapping as following patches

    mirrorY = diag([1,-1, 1, 1]);
    % slot opening iron/air
    nrbRight1 = nrbdegelev(nrbline(point{5}, point{9}), 1);
    nrbRight2 = nrbcirc(params.SR1, point{7}, pi, 1.25*pi);
    nrbRight = nrbglue(nrbRight1, nrbRight2, 2, 1);
    nrbRight.knots = nrbRight.knots/max(nrbRight.knots);
    nrbRight = nrbreverse(nrbRight);
    srf(4) = nrbruled(nrbcirc(params.SRI, center, -beta/2, -gamma1), nrbRight);
    srf(5) = nrbruled(nrbextract(srf(2), 4), nrbline(point{5}, point{6}));
    srf(6) = nrbtform(srf(4), mirrorY);
    % air layer between slot and pin
    nrbRight1 = nrbdegelev(nrbline(point{13}, point{12}), 1);
    nrbRight2 = nrbcirc(params.PR, point{10}, pi, pi+pi/4);
    nrbRight = nrbglue(nrbRight1, nrbRight2, 2, 1);
    nrbRight.knots = nrbRight.knots/max(nrbRight.knots);
    nrbRight = nrbreverse(nrbRight);
    srf(7) = nrbruled(nrbRight, nrbextract(srf(4), 4));
    srf(8) = nrbruled(nrbline(point{13}, point{14}), nrbextract(srf(5), 4));
    srf(9) = nrbtform(srf(7), mirrorY);

    % left part of Pin
    srf(10) = nrbruled(nrbline(point{15}, point{17}), nrbextract(srf(7), 3));
    srf(11) = nrbruled(nrbline(point{17}, point{18}), nrbextract(srf(8), 3));
    srf(12) = nrbruled(nrbline(point{15}, point{16}), nrbextract(srf(11), 3));
    srf(12) = nrbdegelev(srf(12), [0,1]);           % make srf(12) consistent with surrounding patches
    srf(12) = nrbkntins(srf(12), {[], [0.5,0.5]});  % make srf(12) consistent with surrounding patches
    srf(13) = nrbtform(srf(10), mirrorY);
    
    % center bot/top part of pin
    nrbBot1 = nrbcirc(params.SR1, point{7}, 1.25*pi, 1.5*pi);
    nrbBot2 = nrbdegelev(nrbline(point{19}, point{30}), 1);
    nrbBot3 = nrbdegelev(nrbline(point{30}, point{20}), 1);
    nrbBot = nrbglue(nrbglue(nrbBot1, nrbBot2, 2, 1), nrbBot3, 2,1);

    nrbBot.knots = nrbBot.knots/max(nrbBot.knots);

    % nrbBot = nrbreverse(nrbBot);
    nrbTop1 = nrbcirc(params.PR, point{10}, 1.25*pi, 1.5*pi);
    nrbTop2 = nrbdegelev(nrbline(point{21}, point{22}), 1);
    nrbTop12 = nrbglue(nrbTop1, nrbTop2, 2,1);
    nrbTop3 = nrbcirc(params.PR, point{23}, 1.5*pi, 1.75*pi);
    nrbTop = nrbglue(nrbTop12, nrbTop3, 2, 1);
    nrbTop.knots = nrbTop.knots/max(nrbTop.knots);
    srf(14) = nrbruled(nrbBot, nrbTop);
    srf(15) = nrbtform(srf(14), mirrorY);

    srf(16) = nrbruled(nrbline(point{15}, point{25}), nrbextract(srf(14), 4));
    srf(17) = nrbtform(srf(16), mirrorY);
    
    % air between pins
    nrbLeft1 = nrbcirc(params.PR, point{23}, 1.75*pi, 2*pi);
    nrbLeft2 = nrbdegelev(nrbline(point{24}, point{27}), 1);
    nrbLeft3 = nrbcirc(params.PR, point{28}, 0, pi/4);
    nrbLeft = nrbglue(nrbglue(nrbLeft1, nrbLeft2, 2,1), nrbLeft3, 2, 1);
    nrbLeft.knots = nrbLeft.knots/max(nrbLeft.knots);
    nrbRight1 = nrbline(point{20}, point{31});
    nrbRight2 = nrbline(point{31}, point{32});
    nrbRight3 = nrbline(point{32}, point{29});
    nrbRight = nrbglue(nrbglue(nrbRight1, nrbRight2, 2,1), nrbRight3, 2,1);
    nrbRight.knots = nrbRight.knots/max(nrbRight.knots);
    srf(18) = nrbruled(nrbLeft, nrbRight);

    srf(19) = nrbruled(nrbline(point{25}, point{26}), nrbextract(srf(18), 3));
    nrbLeft = nrbextract(srf(12), 3);
    nrbRight = nrbextract(srf(19), 3);
    nrbTop = nrbextract(srf(17), 3);
    nrbBot = nrbextract(srf(16), 3);
    srf(20) = nrbcoons(nrbBot, nrbTop, nrbLeft, nrbRight);
    % iron part above and below first pin
    srf(21) = nrbruled(nrbline(point{4}, point{42}), nrbextract(srf(15), 3));
    srf(22) = nrbtform(srf(21), mirrorY);
    
    % general pin definition
    mirrorX = diag([-1, 1, 1, 1]);
    moveLeft = eye(4);
    moveLeft(1,4) = -point{20}(1);
    moveRight = eye(4);
    moveRight(1,4) = point{5}(1) + params.PP1/2;

    pinSrf(1) = nrbtform(srf(18), moveRight*mirrorX*moveLeft);
    pinSrf(2) = nrbtform(srf(19), moveRight*mirrorX*moveLeft);

    nrbTop1 = nrbline(point{33}, point{34});
    nrbTop2 = nrbline(point{34}, point{35});
    nrbTop3 = nrbline(point{35}, point{29});
    nrbTop = nrbglue(nrbglue(nrbTop1, nrbTop2, 2,1), nrbTop3, 2,1);
    nrbTop.knots = nrbTop.knots/max(nrbTop.knots);

    pinSrf(3) = nrbruled(nrbTop, nrbextract(srf(17), 4)); 
    pinSrf(4) = srf(17);
    pinSrf(5) = nrbtform(pinSrf(3), mirrorY);
    pinSrf(6) = nrbtform(pinSrf(4), mirrorY);
    pinSrf(7) = srf(18);
    pinSrf(8) = srf(19);
    pinSrf(9) = srf(20);
    % copy and move pins
    for iPin = 1:params.PN-2
        nsrf = numel(srf);
        moveRight = eye(4);
        moveRight(1,4) = (params.PH + params.PP1)*iPin;
        for iSrf = 1:numel(pinSrf)
            srf(end+1) = nrbtform(pinSrf(iSrf), moveRight);
        end
        regions.Air = [regions.Air, nsrf + [1,3,5,7]];
        regions.Windings = [regions.Windings, nsrf + [2,4,6,8,9]];
        % patches.Slot(1).Pin(iPin+1).Patches = nsrf + [2,4,6,8,9];
        regions.("Slot1_Pin"+string(iPin+1)) = nsrf + [2,4,6,8,9];
        % Iron parts below and above pin
        line1 = nrbline(point{41+iPin}, point{42+iPin});
        line2 = nrbextract(srf(end-6), 3);
        srf(end+1) = nrbruled(line1, line2);
        srf(end+1) = nrbtform(srf(end), mirrorY);
        regions.Iron = [regions.Iron, nsrf + [10,11]];
    end
    
    % last pin
    nsrf = numel(srf);
    % srf(end+1) = 
    moveRight = eye(4);
    moveRight(1,4) = point{5}(1) + params.PP1/2 + (params.PH + params.PP1)*(params.PN-1);
    
    srf(end+1) = nrbtform(srf(14),moveRight*mirrorX*moveLeft);
    srf(end+1) = nrbtform(srf(16),moveRight*mirrorX*moveLeft);
    srf(end+1) = nrbtform(srf(18),moveRight*mirrorX*moveLeft);
    srf(end+1) = nrbtform(srf(19),moveRight*mirrorX*moveLeft);
    srf(end+1) = nrbtform(srf(15),moveRight*mirrorX*moveLeft);
    srf(end+1) = nrbtform(srf(17),moveRight*mirrorX*moveLeft);
    srf(end+1) = nrbtform(srf(20),moveRight*mirrorX*moveLeft);
    moveRight(1,4) =  (params.PH + params.PP1)*(params.PN-1);
    srf(end+1) = nrbtform(srf(19),moveRight);

    % right side air of last pin
    nrbLeft = nrbextract(srf(19), 4);
    nrbRight1 = nrbcirc(params.SR1, point{40}, 1.75*pi, 2*pi);
    nrbRight2 = nrbdegelev(nrbline(point{38}, point{39}), 1);
    nrbRight3 = nrbcirc(params.SR1, point{41}, 0, pi/4);
    nrbRight = nrbglue(nrbglue(nrbRight1, nrbRight2, 2,1), nrbRight3, 2,1);
    nrbRight.knots = nrbRight.knots/max(nrbRight.knots);
    srfRight = nrbruled(nrbLeft, nrbRight);
    srf(end+1) = nrbtform(srfRight, moveRight);

    regions.Air = [regions.Air, nsrf + [1,3,5,9]];
    regions.Windings = [regions.Windings, nsrf + [2,4,6,7,8]];
    % patches.Slot(1).Pin(params.PN).Patches = nsrf + [2,4,6,7,8];
    regions.("Slot1_Pin"+string(params.PN))= nsrf + [2,4,6,7,8];

    nsrf = numel(srf);
    % last iron part top / bot
    srf(end+1) = nrbruled(nrbline(point{end-1}, point{end}), nrbreverse(nrbextract(srf(end-4), 3)));
    srf(end+1) = nrbtform(srf(end), mirrorY);

    % final iron part right
    srf(end+1) = nrbruled(nrbcirc(params.SRO, [0,0], -alpha, alpha), nrbextract(srf(end-2), 4));
    regions.Iron = [regions.Iron, nsrf + [1,2,3]];
    
    nsrf = numel(srf);
    % Boundaries
    boundaries.Dirichlet.patches = nsrf;
    boundaries.Dirichlet.sides = 3;
    boundaries.Coupling.patches = [1,2,3];
    boundaries.Coupling.sides = [3,3,3];

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
        for iPin = 1:params.PN
            % patches.Slot(I+1).Pin(iPin).Patches = patches.Slot(1).Pin(iPin).Patches + nsrf*I;
            regions.("Slot"+string(I+1)+"_Pin"+string(iPin)) = regions.("Slot1_Pin"+string(iPin)) + nsrf*I;
        end
        boundaries.Dirichlet.patches = [boundaries.Dirichlet.patches, boundariesInit.Dirichlet.patches+nsrf*I];
        boundaries.Dirichlet.sides = [boundaries.Dirichlet.sides, boundariesInit.Dirichlet.sides];
        boundaries.Coupling.patches = [boundaries.Coupling.patches, boundariesInit.Coupling.patches+nsrf*I];
        boundaries.Coupling.sides = [boundaries.Coupling.sides, boundariesInit.Coupling.sides];
    end

    % Rotate again by half a slot and upward facing
    nsrf = numel(srf);
    for I = 1:nsrf
        srf(I) = nrbtform(srf(I), vecrotz(pi/params.SLOTS));
    end

    % Additional Boundaries if not full machine is simulated
    if params.Symmetry ~= 1
        boundaries.Periodic1.patches = [1, 11+11*(1:params.PN)];
        boundaries.Periodic1.sides = [1,repmat(3, 1, params.PN)];
        boundaries.Periodic2.patches = nsrf - (params.PN+1)*11 + [2, 9+11*(1:params.PN)];
        boundaries.Periodic2.sides = [2, repmat(3, 1, params.PN)];
    end
    if params.LumpWindings
        regionsNoSlots.Iron = regions.Iron;
        regionsNoSlots.Air = regions.Air;
        regionsNoSlots.Windings = regions.Windings;
        regions = regionsNoSlots;
    else
        regions = rmfield(regions, "Windings");
    end

    if (params.preview)
        preview_nurbs(srf, regions);
        % for i = 1:numel(point)
        %     scatter(point{i}(1), point{i}(2), "black", "filled");
        %     text(point{i}(1)+1e-5, point{i}(2), string(i))
        % end
        view(2)
        axis equal
    end
end
   