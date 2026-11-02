% https://www.mathcha.io/editor/Z0NnmfdrHpVHKmnWJ6H7Mzvejs7gxVEViw1nJGM

function [srf, regions, boundaries, params] = GEO_DISstator_1(params)
    arguments
        params.preview (1,1) logical = true;
        params.Symmetry = 4;
        params.SLOTS = 36;             % Nuber of stator slots
        
        params.HSO = 1.1e-4;           % Height Slot Opening
        params.HY = 20e-3;             % Yoke Height
        params.RAG = 62.75e-3;         % Radius Air Gap
        params.SRI = 63e-3;            % Stator Radius Inner
        params.SRO = 100e-3;           % Stator Radius Outer
        params.SR1 = 1.2e-3;           % Slot Radius
        params.WSO = 2.5e-3;           % Width Slot Opening
        params.WT = 6e-3;              % Width Tooth
        params.LumpWindings (1,1) logical = false;
    end
    
    % calculated parameters
    center = [0, 0, 0];    
    beta = 2*pi/params.SLOTS;
    alpha = beta/2;
    gamma1 = 2*asin(params.WSO/2 / params.SRI);
    L2 = params.SRO - params.HY - params.SR1;           % lenght till center of circle 3
    
    % Equations:
    % L1 = cos(gamma)*params.SRI + params.SW3 + Dslit
    % (params.SW2/2)^2 + Dslit^2 = Sr2^2
    % sin(alpha)*L1 = params.SW1/2 + Sr2
    const = params.WT/2/sin(alpha) - cos(gamma1/2)*params.SRI - params.HSO;
    a = (1/sin(alpha)^2-1);
    b = (2/sin(alpha)*const);
    c = const^2 + (params.WSO/2)^2;
    Sr2 = (-b + sqrt(b^2-4*a*c))/2/a;
    Dslot = (Sr2^2 - (params.WSO/2)^2)^0.5;
    L1 = (params.WT/2 + Sr2)/sin(alpha);
    
    Lslot = sqrt(Sr2^2 - params.WSO^2/4);
    Rslot = sqrt((L1-Lslot)^2 + (params.WSO/2)^2);
    gamma2 = 2*asin(params.WSO/2 / Rslot);
    epsilon = asin(params.WSO/2/Sr2);
    
    B = params.WT/2 + params.SR1;
    angleL2 = atan((L2*sin(alpha)-B)/(L2*cos(alpha)));
    L3 = L2*tan(angleL2);
    
    %setting points
    point{1} = [params.SRI*cos(alpha), -params.SRI*sin(alpha)];
    point{2} = [params.SRI*cos(gamma1/2), -params.SRI*sin(gamma1/2)];
    point{3} = [params.SRI*cos(gamma1/2), params.SRI*sin(gamma1/2)];
    point{4} = [params.SRI*cos(alpha), params.SRI*sin(alpha)];
    
    point{5} = [Rslot*cos(alpha), -Rslot*sin(alpha)];
    point{6} = [L1-Lslot, - params.WSO/2];
    point{7} = [L1-Lslot, params.WSO/2];
    point{8} = [Rslot*cos(alpha), Rslot*sin(alpha)];
    
    point{9} = [L1 - sin(alpha)*Sr2, -cos(alpha)*Sr2];
    point{10} = [L1, 0];
    point{11} = [L1 - sin(alpha)*Sr2, cos(alpha)*Sr2];
    
    point{12} = [L2 - sin(alpha)*params.SR1, -L3 - cos(alpha)*params.SR1];
    point{13} = [L2, -L3];
    point{14} = [L2, L3];
    point{15} = [L2 - sin(alpha)*params.SR1, L3 + cos(alpha)*params.SR1];
    
    point{16} = [L2+params.SR1, -L3];
    point{17} = [L2+params.SR1, L3];
    
    point{18} = params.SRO*[cos(alpha), -sin(alpha)];
    point{19} = params.SRO*[cos(alpha), sin(alpha)];
    
    L4 =  (point{17}(1) - point{10}(1));
    beta2 = pi/2-alpha-epsilon;
    beta4 = pi/2 + alpha;
    a1 = 0.5*Dslot*params.WSO/2;
    a2 = Sr2^2/2*beta2;
    a3 = 1/2*Sr2^2*tan(alpha);
    As = params.SR1^2*tan(beta4/2) - beta4/2*params.SR1^2;
    % Solve quadratic formula for having same slot areas:
    a = tan(alpha);
    b = 2*Sr2/cos(alpha);
    c = a1+a2+a3+As - L4*Sr2/cos(alpha) - 0.5*L4^2*tan(alpha);
    Ll = (-b + (b^2- 4*a*c)^0.5)/(2*a);
    Lr = L4 - Ll;
    % Double check areas:
    a4 = Sr2/cos(alpha)*Ll;
    a5 = 1/2*Ll^2*tan(alpha);
    Aleft = a1+a2+a3+a4+a5;
    Aright = Lr*(Sr2/cos(alpha)+tan(alpha)*Ll) + 1/2*Lr*Lr*tan(alpha) - As;
    
    point{21} = [point{10}(1) + Ll, -Sr2/cos(alpha) - Ll*tan(alpha)]; % set for equal slot size
    point{22} = [point{10}(1) + Ll, Sr2/cos(alpha) + Ll*tan(alpha)]; % set for equal slot size
    
    point{20} = (point{5}+point{18})/2;
    point{23} = (point{8}+point{19})/2;
    point{20} = point{20} * vecnorm(point{21})/vecnorm(point{20});
    point{23} = point{23} * vecnorm(point{22})/vecnorm(point{23});
    
    pS = 0.6; % point share > 0.5 (centered)
    point{24} = pS*pS*point{9} + pS*(1-pS)*(point{11}) + pS*(1-pS)*point{21} + (1-pS)*(1-pS)*point{22};
    point{25} = pS*pS*point{11} + pS*(1-pS)*(point{9}) + pS*(1-pS)*point{22} + (1-pS)*(1-pS)*point{21};
    
    point{26} = pS*pS*point{12} + pS*(1-pS)*(point{15}) + pS*(1-pS)*point{21} + (1-pS)*(1-pS)*point{22};
    point{27} = pS*pS*point{15} + pS*(1-pS)*(point{22}) + pS*(1-pS)*point{12} + (1-pS)*(1-pS)*point{21};
    
    point{28} = pS*point{21}+ (1-pS)*point{22};
    point{29} = pS*point{22}+ (1-pS)*point{21};
    
    % parameter checks for motor:
    L1test = point{10}(1)*1000;
    R2test = Sr2*1000;
    L2test = point{13}(1)*1000;
    Wtest = point{14}(2)*2000;
    
    % creating surfaces
    % Air inside COUPLING INTERFACE IS HARD CODED FOR SAME BOUNDARY REPRESENTATION
    srf(1) = nrbruled(nrbcirc(params.RAG, center, -beta/2, -beta/6), nrbcirc(params.SRI, center, -beta/2, -gamma1/2));
    srf(2) = nrbruled(nrbcirc(params.RAG, center, -beta/6, beta/6), nrbcirc(params.SRI, center, -gamma1/2, gamma1/2));
    srf(3) = nrbruled(nrbcirc(params.RAG, center, beta/6, beta/2), nrbcirc(params.SRI, center, gamma1/2, beta/2));
    % Iron/air layer
    srf(4) = nrbruled(nrbcirc(params.SRI, center, -beta/2, -gamma1/2), nrbcirc(Rslot, center, -beta/2, -gamma2/2));
    srf(5) = nrbruled(nrbcirc(params.SRI, center, -gamma1/2, gamma1/2), nrbline(point{6}, point{7}));
    srf(6) = nrbruled(nrbcirc(params.SRI, center, gamma1/2, beta/2), nrbcirc(Rslot, center, gamma2/2, beta/2));
    
    % Copper Slot Top Left
    nrbTop1 = nrbdegelev(nrbline(point{22}, point{11}), 1);
    nrbTop2 = nrbcirc(Sr2, point{10}, pi/2+alpha, pi-epsilon);
    nrbTop = nrbglue(nrbTop1, nrbTop2, 2, 1);
    nrbTop.knots = nrbTop.knots/max(nrbTop.knots);
    nrbTop = nrbreverse(nrbTop);
    nrbLeft = nrbline(point{25}, point{7});
    nrbRight = nrbline(point{29}, point{22});
    nrbBot = nrbline(point{25}, point{29});
    srf(7) = nrbcoons(nrbBot, nrbTop, nrbLeft, nrbRight);
    
    % Copper Slot Top Right
    nrbTop1 = nrbdegelev(nrbline(point{15}, point{22}), 1);
    nrbTop2 = nrbcirc(params.SR1, point{14}, 0, pi/2+alpha);
    nrbTop = nrbglue(nrbTop1, nrbTop2, 1, 2);
    nrbTop.knots = nrbTop.knots/max(nrbTop.knots);
    nrbTop = nrbreverse(nrbTop);
    nrbLeft = nrbextract(srf(7), 2);
    nrbBot = nrbline(point{29}, point{27});
    nrbRight = nrbline(point{27}, point{17});
    srf(8) = nrbcoons(nrbBot, nrbTop, nrbLeft, nrbRight);
    
    % Copper Slot Bot Left and Right
    mirror = eye(4);
    mirror(2, 2) = -1;
    srf(9) = nrbtform(srf(7), mirror);
    srf(10) = nrbtform(srf(8), mirror);
    
    % Copper Slot Left
    nrbBot = nrbreverse(nrbextract(srf(9), 1));
    nrbTop = nrbreverse(nrbextract(srf(7), 1));
    nrbLeft = nrbextract(srf(5), 4);
    nrbRight = nrbline(point{24}, point{25});
    srf(11) = nrbcoons(nrbBot, nrbTop, nrbLeft, nrbRight);
    %
    % Copper Slot Center Left
    nrbBot = nrbextract(srf(9), 3);
    nrbTop = nrbextract(srf(7), 3);
    nrbLeft = nrbextract(srf(11), 2);
    nrbRight = nrbline(point{28}, point{29});
    srf(12) = nrbcoons(nrbBot, nrbTop, nrbLeft, nrbRight);
    
    % Copper Slot Center Right
    nrbBot = nrbextract(srf(10), 3);
    nrbTop = nrbextract(srf(8), 3);
    nrbLeft = nrbextract(srf(12), 2);
    nrbRight = nrbline(point{26}, point{27});
    srf(13) = nrbcoons(nrbBot, nrbTop, nrbLeft, nrbRight);
    
    % Copper Slot Right
    nrbBot = nrbextract(srf(10), 2);
    nrbTop = nrbextract(srf(8), 2);
    nrbLeft = nrbextract(srf(13), 2);
    nrbRight = nrbline(point{16}, point{17});
    srf(14) = nrbcoons(nrbBot, nrbTop, nrbLeft, nrbRight);
    
    % Iron Top Left
    nrbBot = nrbextract(srf(7), 4);
    nrbTop = nrbline(point{8}, point{23});
    nrbLeft = nrbextract(srf(6), 4);
    nrbRight = nrbline(point{22}, point{23});
    srf(15) = nrbcoons(nrbBot, nrbTop, nrbLeft, nrbRight);
    
    % Iron Top Right
    nrbBot = nrbextract(srf(8), 4);
    nrbTop = nrbline(point{23}, point{19});
    nrbLeft = nrbextract(srf(15), 2);
    nrbRight = nrbline(point{17}, point{19});
    srf(16) = nrbcoons(nrbBot, nrbTop, nrbLeft, nrbRight);
    
    % Bot iron
    srf(17) = nrbtform(srf(15), mirror);
    srf(18) = nrbtform(srf(16), mirror);
    
    % Right iron
    nrbBot = nrbextract(srf(18), 2);
    nrbTop = nrbextract(srf(16), 2);
    nrbLeft = nrbextract(srf(14), 2);
    nrbRight = nrbcirc(params.SRO, center, -alpha, alpha);
    srf(19) = nrbcoons(nrbBot, nrbTop, nrbLeft, nrbRight);

    % Rotate again by half a slot and upward facing
    nsrf = numel(srf);
    for I = 1:nsrf
        srf(I) = nrbtform(srf(I), vecrotz(pi/params.SLOTS));
    end

    % regions definitions
    regions.Iron = [4, 6, 15:19];
    regions.Air = [1:3, 5];
    regions.Slot1up = [8,10,13,14];
    regions.Slot1down = [7,9,11,12];
    regions.Windings = [regions.Slot1up, regions.Slot1down];

    % boundaries
    boundaries.Dirichlet.patches = 19;
    boundaries.Dirichlet.sides = 2;
    boundaries.Coupling.patches = [1,2,3];
    boundaries.Coupling.sides = [3,3,3];
    
    % Copy and rotate
    nsrf = numel(srf);
    IronInit = regions.Iron;
    AirInit = regions.Air;
    boundariesInit = boundaries;
    WindingsInit = regions.Windings;
    for I = 1:params.SLOTS/params.Symmetry-1
        for i = 1:nsrf
            srf(I*nsrf+i) = nrbtform(srf(i), vecrotz(I*2*pi/params.SLOTS));
        end
        regions.Iron = [regions.Iron, IronInit+nsrf*I];
        regions.Air = [regions.Air, AirInit+nsrf*I];
        regions.("Slot"+num2str(I+1)+"up") = regions.Slot1up + nsrf*I;
        regions.("Slot"+num2str(I+1)+"down") = regions.Slot1down + nsrf*I;
        regions.Windings = [regions.Windings, WindingsInit + nsrf*I];
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
    if params.Symmetry ~= 1
        boundaries.Periodic1.patches = [1,4,17,18];
        boundaries.Periodic1.sides = [1,1,4,4];
        boundaries.Periodic2.patches = nsrf-[16,13,4,3];
        boundaries.Periodic2.sides = [2,2,4,4];
    end
    
    if (params.preview)
        preview_nurbs(srf, regions);
        % for i = 1:numel(point)
        %     scatter(point{i}(1), point{i}(2), "black", "filled");
        %     text(point{i}(1)+1e-5, point{i}(2), string(i))
        % end
    end
end
   
