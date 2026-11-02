

function c = GEO_IPMrotor_V_1_constr(params)
    arguments
        params.preview = false;
        params.POLES = 4;
        params.Symmetry = 4;
        params.REFINEMENT = [1, 0, 3, 0, 1];
    
        params.RRO = 50e-3;             % Rotor Radius Outer
        params.RRI = 14.55e-3;          % Rotor Radius Inner
        params.RAG = 50.5e-3;           % Radius Air Gap
    
        params.DMAG = 30e-3;            % Distance Magnet Positioning 
        params.LSLIT1 = 6.4e-3;         % Air Slit Positioning
        params.LSLIT2 = 4.3e-3;         % Air Slit Positioning
        params.DSLIT5 = 1.0e-3;         % Air Slit Positioning
        params.DSLIT6 = 2.0e-3;         % Air Slit Positioning
        params.MA = 150;                % Magnet opening Angle
        params.MH = 4.0e-3;             % Magnet Height
        params.MW = 22.0e-3;            % Magnet Width
        params.RA1 = 144.2;             % Air Slit Angle 1
        params.RA2 = 166.0;             % Air Slit Angle 2
        params.RF = 0.2e-3;             % Radius Fillets 
        params.RS = 1.0e-3;             % Rotor Iron Slit thickness
        params.RW2 = 1.0e-3;            % Air Slit Positioning
        params.RW3 = 1.0e-3;            % Air Slit Positioning
        params.RW4 = 1.0e-3;            % Air Slit Positioning
        params.RW5 = 1.0e-3;            % Air Slit Positioning
        params.WMAG = 4.0e-3;           % Width Magnet Positoning
        
        params.LumpMagnets (1,1) logical = false;
        params.minIron = 1.5e-3;        % Constraints for iron bridges
    end

    % calculated parameters
    beta = pi/params.POLES;
    gamma = atan(params.WMAG/params.DMAG);
    delta = beta-gamma;
    epsilon = deg2rad(params.MA/2)-beta;
    phi1 = deg2rad(rad2deg(epsilon) + params.RA1 - 180);
    phi2 = deg2rad(rad2deg(epsilon) + params.RA2 - 180);

    % point definitions
    point{1} = [cos(beta), sin(beta)]*params.RRI;
    % magnet
    point{2} = [cos(delta), sin(delta)]*(params.DMAG^2+params.WMAG^2)^0.5;
    point{3} = point{2} + [cos(epsilon), -sin(epsilon)]*params.MW;
    point{4} = point{2} + [-sin(epsilon), -cos(epsilon)]*params.MH;
    point{5} = point{3} + [-sin(epsilon), -cos(epsilon)]*params.MH;
    % rw points
    point{6} = point{2} + [-cos(epsilon), sin(epsilon)]*params.RW4; 
    point{7} = point{3} + [cos(epsilon), -sin(epsilon)]*params.RW2; 
    point{8} = point{4} + [-cos(epsilon), sin(epsilon)]*params.RW5; 
    point{9} = point{5} + [cos(epsilon), -sin(epsilon)]*params.RW3; 

    point{10} = [(params.DMAG-params.DSLIT5)*cos(beta), (params.DMAG-params.DSLIT5)*sin(beta)];
    point{11} = [(params.DMAG-params.DSLIT6)*cos(beta), (params.DMAG-params.DSLIT6)*sin(beta)];
    
    point{12} = point{10} + [sin(beta), -cos(beta)]*params.RS/2;
    point{13} = point{11} + [sin(beta), -cos(beta)]*params.RS/2;

    point{14} = point{9} + [cos(phi1), -sin(phi1)]*params.LSLIT1;
    point{15} = point{7} + [cos(phi2), -sin(phi2)]*params.LSLIT2;

    point{16} = [params.RRI, 0];
    point{17} = [point{5}(1), 0];
    point{18} = [params.RRO, 0];
    
    point{19} = vecnormalize(point{15})*params.RRO;
    point{20} = vecnormalize(point{3})*params.RRO;
    point{21} = vecnormalize(point{2})*params.RRO;
    point{22} = vecnormalize(point{13})*params.RRO;

    point{23} = point{10} - [sin(beta), -cos(beta)]*params.RS/2;
    point{24} = point{11} - [sin(beta), -cos(beta)]*params.RS/2;
    point{25} = vecnormalize(point{13})*params.RRI;

    c(1) = vecmag(point{9}) - (params.RRO) + params.minIron;
    c(2) = vecmag(point{14}) - (params.RRO) + params.minIron;
    c(3) = -point{14}(2) + params.minIron/2;
    c(4) = -point{9}(2) + params.minIron/2;
    c(5) = vecmag(point{15}) - (params.RRO) + params.minIron;
    c(6) = -vecmag(point{15}-point{14}) + 2*params.RF;
    c(7) = -vecmag(point{12}-point{13}) + 2*params.RF;
    c(8) = params.RRI + 2*params.minIron -vecmag(point{8}) ;
    c(9) = params.RRI + 2*params.minIron -vecmag(point{13});

    c = c * 100;    % scale to create better convergence behavior
end

