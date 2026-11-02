% https://www.mathcha.io/editor/BY534Czkt6JsJ8Fx9BY0gCVrPq8XtVplvWoFJEmk54

function [srf, regions, boundaries, params] = GEO_IMrotor_C_1(params)
    arguments
        params.preview (1,1) logical = true;
        params.Symmetry = 4;
        params.BARS = 28;              % Number of Rotor Bars
        
        params.RAG = 62.75e-3;         % Radius Air Gap
        params.RRI = 20e-3;            % Rotor Radius Inner
        params.RRO = 62.5e-3;          % Rotor Radius Outer
        
        params.RR1 = 41.49e-3;          % Rotor Radius 1
        params.RR2 = 51.48e-3;          % Rotor Radius 2
        params.RR3 = 58.95e-3;          % Rotor Radius 3
        params.RR4 = 0.99e-3;           % Rotor Radius 4
        params.RR5 = 2.11e-3;           % Rotor Radius 5
        params.RR6 = 5.7e-3/2;          % Rotor Radius 6
        params.RWSO = 1e-3;             % Rotor Width Slot Opening
        params.LumpBars (1,1) logical = false;
    end
    
    center = [0, 0, 0];

    % patch definitions
    regions.Iron = [1:3, 9, 10, 12,13,15,16];
    regions.Air = [18,19,20];
    regions.Aluminum = [4:8,11,14,17];
    regions.Bar1 = regions.Aluminum;

    % calculated parameters
    beta = 2*pi/params.BARS;
    alpha = beta/2;
    gamma = acos((params.RR5-params.RR4)/(params.RR2-params.RR1));
    delta = pi - gamma; 

    ank7 = params.RR1 + params.RR4*cos(delta);
    gek7 = params.RR4*sin(delta);
    rad7 = sqrt(ank7^2 + gek7^2);
    ang7 = atan(gek7/ank7); % winkel an punkt 7

    ank17 = params.RR2 + params.RR5*cos(delta);
    gek17 = params.RR5*sin(delta);
    rad17 = sqrt(ank17^2 + gek17^2);
    ang17 = atan(gek17/ank17); % winkel an punkt 17

    ank15 = sqrt(params.RRO^2 - (params.RWSO/2)^2);
    ang15 = atan((params.RWSO/2)/ank15); % winkel an punkt 15
    ang18 = atan((params.RWSO/2)/(sqrt(params.RR5^2 - (params.RWSO/2)^2))); % winkel an 18
    
    rm2 = params.RR3 + sqrt(params.RR6^2 - (params.RWSO/2)^2);
    rm2_rot = sqrt(rm2^2 + (params.RWSO/2)^2);  %radius an punkt 11
    ang11 = asin((params.RWSO/2)/rm2_rot);   %winkel an punkt 11

    %setting points
    point{1} = [params.RRI*cos(alpha), -params.RRI*sin(alpha)];
    point{2} = [params.RRI*cos(1/3*alpha), -params.RRI*sin(1/3*alpha)];
    point{3} = [params.RRI*cos(1/3*alpha), params.RRI*sin(1/3*alpha)];
    point{4} = [params.RRI*cos(alpha), params.RRI*sin(alpha)];

    point{5} = [rad7*cos(alpha), -rad7*sin(alpha)];
    point{6} = [rad7*cos(ang7), -rad7*sin(ang7)];
    point{7} = [rad7*cos(ang7), rad7*sin(ang7)];
    point{8} = [rad7*cos(alpha), rad7*sin(alpha)];

    point{9} = [(rm2_rot)*cos(alpha), -(rm2_rot)*sin(alpha)];
    point{10} = [params.RR3 + sqrt(params.RR6^2-(params.RWSO/2)^2), -(params.RWSO/2)];
    point{11} = [params.RR3 + sqrt(params.RR6^2-(params.RWSO/2)^2), (params.RWSO/2)];
    point{12} = [(rm2_rot)*cos(alpha), (rm2_rot)*sin(alpha)];

    point{13} = [params.RRO*cos(alpha), -params.RRO*sin(alpha)];
    point{14} = [sqrt(params.RRO^2 - (params.RWSO/2)^2), -(params.RWSO/2)];
    point{15} = [sqrt(params.RRO^2 - (params.RWSO/2)^2), (params.RWSO/2)];
    point{16} = [params.RRO*cos(alpha), params.RRO*sin(alpha)];

    point{17} = [rad17*cos(ang17), rad17*sin(ang17)];
    point{18} = [params.RR2 + params.RR5*cos(ang18), (params.RWSO/2)];
    point{19} = [params.RR2 + params.RR5*cos(ang18), -(params.RWSO/2)];
        
    point{20} = vecnorm(point{18}) * [cos(alpha), sin(alpha)];

    point{21} = [point{7}(1)*0.9 + point{18}(1)*0.1, point{7}(2)*0.6 + point{6}(2)*0.4];
    point{22} = [point{6}(1)*0.9 + point{19}(1)*0.1, point{7}(2)*0.4 + point{6}(2)*0.6];

    point{23} = [point{18}(1)*0.9 + point{7}(1)*0.1, point{18}(2)*0.6 + point{19}(2)*0.4];
    point{24} = [point{18}(1)*0.9 + point{7}(1)*0.1, point{18}(2)*0.4 + point{19}(2)*0.6];

    point{25} = [params.RR3 - sqrt(params.RR6^2-(params.RWSO/2)^2), (params.RWSO/2)];
    point{26} = [params.RR3 - sqrt(params.RR6^2-(params.RWSO/2)^2), -(params.RWSO/2)];

    point{27} = vecnorm(point{25}) * [cos(alpha), sin(alpha)];
    point{28} = vecnorm(point{25}) * [cos(alpha), -sin(alpha)];

    % create nurbs
    % Top Iron
    line1 = nrbdegelev(nrbline(point{7}, point{17}), 1);
    line2 = nrbcirc(params.RR5, [params.RR2, 0, 0], ang18, delta);
    nrbBot = nrbglue(line1, line2);
    nrbBot.knots = nrbBot.knots / max(nrbBot.knots);
    nrbLeft = nrbline(point{4}, point{7});
    nrbTop = nrbline(point{4}, point{20});
    nrbRight = nrbcirc(vecnorm(point{18}), center, atan2(point{18}(2), point{18}(1)), alpha);
    helpsrf = nrbcoons(nrbBot, nrbTop, nrbreverse(nrbLeft), nrbRight);
    % manual corrections of control points
    p1 = helpsrf.coefs(1:2, 1, end);    % space points on lower boundary correctly
    p2 = helpsrf.coefs(1:2, end, end);
%     pos = [0.5, 0.25, 0.15, 0.1, 0.075, 0.05, 0.025];
    pos = [0.4, 0.2, 0.1];
    helpsrf.coefs(1:2, 2:end-1, end) = p1.*pos + p2.*(1-pos);
    srf(1) = correctInnerPoints(helpsrf, 'x');

    % Bottom Iron
    mirror = eye(4);
    mirror(2, 2) = -1;
    srf(2) = nrbtform(srf(1), mirror);

    % Left Iron
    leftCirc = nrbcirc(params.RRI, center, -alpha, alpha);
    rightCirc = nrbreverse(nrbcirc(params.RR4, [params.RR1, 0, 0], delta, 2*pi-delta));
    nrbTop = nrbreverse(nrbextract(srf(1), 1));
    nrbBot = nrbreverse(nrbextract(srf(2), 1));
    srf(3) = nrbcoons(nrbBot, nrbTop , leftCirc, rightCirc);
    
    % Left Copper
    rightSide = nrbline(point{22}, point{21});
    leftSide = nrbextract(srf(3), 2);
    srf(4) = nrbruled(leftSide, rightSide);

    % Top Copper
    LeftSide = nrbreverse(nrbextract(srf(4), 2));
    TopSide = nrbextract(srf(1), 3);
    RightSide = nrbline(point{23} ,point{18});
    BotSide = nrbline(point{21} ,point{23});
    srf(5) = nrbcoons(BotSide, TopSide, LeftSide, RightSide);
    % Bot Copper
    srf(6) = nrbtform(srf(5), mirror);
    % Center Copper
    LeftSide = nrbextract(srf(4), 4);
    TopSide = nrbextract(srf(5), 3);
    RightSide = nrbline(point{24} ,point{23});
    BotSide = nrbextract(srf(6), 3);
    srf(7) = nrbcoons(BotSide, TopSide, LeftSide, RightSide);
    
    % Center right copper
    LeftSide = nrbextract(srf(7), 2);
    TopSide = nrbextract(srf(5), 2);
    RightSide = nrbline(point{19} ,point{18});
    BotSide = nrbextract(srf(6), 2);
    srf(8) = nrbcoons(BotSide, TopSide, LeftSide, RightSide);

    % Top and center Iron
    RightSide = nrbcirc(vecnorm(point{25}), center, atan2(point{25}(2), point{25}(1)), alpha);
    LeftSide = nrbextract(srf(1), 2);
    srf(9) = nrbruled(LeftSide, RightSide);
    srf(10) = nrbtform(srf(9), mirror);

    % Center copper
    LeftSide = nrbextract(srf(8), 2);
    TopSide = nrbextract(srf(9), 1);
    BotSide = nrbextract(srf(10), 1);
    RightSide = nrbline(point{26}, point{25});
    srf(11) = nrbcoons(BotSide, TopSide, LeftSide, RightSide);

    % Top and Bot Iron 
    LeftSide = nrbextract(srf(9), 4);
    RightSide = nrbcirc(vecnorm(point{11}), center, atan2(point{11}(2), point{11}(1)), alpha);
    TopSide = nrbline(point{27}, point{12});
    BotSide = nrbreverse(nrbcirc(params.RR6, [params.RR3, 0, 0], asin((params.RWSO/2)/params.RR6), pi - asin((params.RWSO/2)/params.RR6)));
    srf(12) = nrbcoons(BotSide, TopSide, LeftSide, RightSide);
    srf(13) = nrbtform(srf(12), mirror);

    % Round center Copper
    LeftSide = nrbextract(srf(11), 2);
    RightSide = nrbline(point{10}, point{11});
    TopSide = nrbextract(srf(12), 3);
    BotSide = nrbextract(srf(13), 3);
    srf(14) = nrbcoons(BotSide, TopSide, LeftSide, RightSide);

    % Very right iron
    RightSide = nrbcirc(vecnorm(point{15}), center, atan2(point{15}(2), point{15}(1)), alpha);
    LeftSide = nrbextract(srf(12), 2);
    srf(15) = nrbruled(LeftSide, RightSide);
    srf(16) = nrbtform(srf(15), mirror);
    
    % Very right copper
    LeftSide = nrbextract(srf(14), 2);
    RightSide = nrbcirc(vecnorm(point{15}), center, atan2(point{14}(2), point{14}(1)), atan2(point{15}(2), point{15}(1)));
    TopSide = nrbextract(srf(15), 1);
    BotSide = nrbextract(srf(16), 1);
    srf(17) = nrbcoons(BotSide, TopSide, LeftSide, RightSide);

    % Air Slit top and bot
    centerAngle = beta/4;
    LeftSide = nrbextract(srf(15), 4);
    RightSide = nrbcirc(params.RAG, center, centerAngle/2, alpha);
    srf(18) = nrbruled(LeftSide, RightSide);
    srf(19) = nrbtform(srf(18), mirror);

    % Air slit center
    LeftSide = nrbextract(srf(17), 2);
    RightSide = nrbcirc(params.RAG, center, -centerAngle/2, centerAngle/2);
    srf(20) = nrbruled(LeftSide, RightSide);

    % Boundaries
    boundaries.Dirichlet.patches = 3;
    boundaries.Dirichlet.sides = 1;
    boundaries.Coupling.patches = [18,19,20];
    boundaries.Coupling.sides = [4,4,4];

    % Copy and rotate
    nsrf = numel(srf);
    IronInit = regions.Iron;
    AirInit = regions.Air;
    boundariesInit = boundaries;
    for I = 1:params.BARS/params.Symmetry-1
        for i = 1:nsrf
            srf(I*nsrf+i) = nrbtform(srf(i), vecrotz(I*pi/params.BARS*2));
        end
        regions.Iron = [regions.Iron, IronInit+nsrf*I];
        regions.Air = [regions.Air, AirInit+nsrf*I];
        regions.("Bar"+num2str(I+1)) = regions.Bar1+nsrf*I;
        regions.Aluminum = [regions.Aluminum, regions.("Bar"+num2str(I+1))];
        boundaries.Dirichlet.patches = [boundaries.Dirichlet.patches, boundariesInit.Dirichlet.patches+nsrf*I];
        boundaries.Dirichlet.sides = [boundaries.Dirichlet.sides, boundariesInit.Dirichlet.sides];
        boundaries.Coupling.patches = [boundaries.Coupling.patches, boundariesInit.Coupling.patches+nsrf*I];
        boundaries.Coupling.sides = [boundaries.Coupling.sides, boundariesInit.Coupling.sides];
    end
    if params.LumpBars
        regionsNoBars.Iron = regions.Iron;
        regionsNoBars.Air = regions.Air;
        regionsNoBars.Aluminum = regions.Aluminum;
        regions = regionsNoBars;
    else
        regions = rmfield(regions, "Aluminum");
    end

    nsrf = numel(srf);
    for I = 1:nsrf
        srf(I) = nrbtform(srf(I), vecrotz(pi/params.BARS));
    end

    % Additional Boundaries if not full machine is simulated
    if params.Symmetry == 1
        % do nothing
    else
        boundaries.Periodic1.patches = [2,10,13,16,19];
        boundaries.Periodic1.sides = [4,2,4,2,2];
        boundaries.Periodic2.patches = nsrf-[19,11,8,5,2];
        boundaries.Periodic2.sides = [4,2,4,2,2];
    end

    if (params.preview)
        preview_nurbs(srf, regions);
        % for i = 1:numel(point)
        %     scatter(point{i}(1), point{i}(2), "black", "filled");
        %     text(point{i}(1)+1e-5, point{i}(2), string(i))
        % end
    end
end

function srf = correctInnerPoints(srf, dim)
    if dim == 'x'
        for i = 2:size(srf.coefs, 2)-1
            x = linspace(srf.coefs(1, i, 1)./srf.coefs(4, i, 1), srf.coefs(1, i, end)./srf.coefs(4, i, end), size(srf.coefs, 3));
            y = linspace(srf.coefs(2, i, 1)./srf.coefs(4, i, 1), srf.coefs(2, i, end)./srf.coefs(4, i, end), size(srf.coefs, 3));
            srf.coefs(1, i, :) = squeeze(x).*squeeze(srf.coefs(4, i, :))';
            srf.coefs(2, i, :) = squeeze(y).*squeeze(srf.coefs(4, i, :))';
        end
    elseif dim == 'y'
        for i = 2:size(srf.coefs, 3)-1
            x = linspace(srf.coefs(1, 1, i)./srf.coefs(4, 1, i), srf.coefs(1, end, i)./srf.coefs(4, end, i), size(srf.coefs, 2));
            y = linspace(srf.coefs(2, 1, i)./srf.coefs(4, 1, i), srf.coefs(2, end, i)./srf.coefs(4, end, i), size(srf.coefs, 2));
            srf.coefs(1, :, i) = squeeze(x).*squeeze(srf.coefs(4, :, i));
            srf.coefs(2, :, i) = squeeze(y).*squeeze(srf.coefs(4, :, i));
        end
    end
end
