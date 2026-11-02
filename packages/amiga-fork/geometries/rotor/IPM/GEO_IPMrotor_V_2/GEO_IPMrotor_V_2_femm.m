

function [params, magnetization] = GEO_IPMrotor_V_2_femm (params)
    arguments
        params.preview = true;
        params.POLES = 8;
        params.Symmetry = 8;
        params.REFINEMENT = [0, 0, 3, 1];
        params.MeshSize = 0;      % Automesh
        params.MeshSizeAG = 1e-3;    % mesh sides at air gap
    
        params.RRO = 78.6e-3;             % Rotor Radius Outer
        params.RRI = 26.5e-3;          % Rotor Radius Inner
        params.RAG = 78.9e-3;           % Radius Air Gap
    
        params.DMAG = 70.9e-3;            % Distance Magnet Positioning 
        params.LSLIT1 = 1.9e-3;         % Air Slit Positioning
        params.LSLIT2 = 6e-3;           % Air Slit Positioning
        params.DSLIT1 = 3.5e-3;         % Air Slit Positioning
        params.DSLIT2 = 0.1e-3;         % Air Slit Positioning
        params.MA = 160;                % Magnet opening Angle
        params.MH = 4.0e-3;             % Magnet Height
        params.MW = 22.0e-3;            % Magnet Width
        params.RA1 = 190;               % Air Slit Angle 1
        params.RA2 = 210.0;             % Air Slit Angle 2
        params.RS = 1.4e-3;             % Rotor Iron Slit thickness
        params.WMAG = 1.7e-3;           % Width Magnet Positoning
    end
    
    mytimer("FEMM", "start");

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

    point{12} = params.RRI*[cos(motorAngle/2), -sin(motorAngle/2)];
    point{13} = params.RRO*[cos(motorAngle/2), -sin(motorAngle/2)];
    point{14} = params.RAG*[cos(motorAngle/2), -sin(motorAngle/2)];
    point{15} = params.RRI*[cos(motorAngle/2), sin(motorAngle/2)];
    point{16} = params.RRO*[cos(motorAngle/2), sin(motorAngle/2)];
    point{17} = params.RAG*[cos(motorAngle/2), sin(motorAngle/2)];

    % mirror upper magnet
    for iPoint = 1:10
        point{17+iPoint} = [point{iPoint}(1), -point{iPoint}(2)];
    end

    % Rotate all points
    for iPoint = 1:numel(point)
        point{iPoint} = [cos(motorAngle/2), -sin(motorAngle/2); sin(motorAngle/2), cos(motorAngle/2)]*point{iPoint}';
    end

    openfemm();
    newdocument(0);
    mi_probdef(0, 'meters', 'planar' ,1e-6, 0.1 ,30, 0);

    % GENERATE GEOMETRY
    % outer line
    mi_drawline(point{12}, point{13});
    mi_drawline(point{13}, point{14});
    mi_drawarc(point{14}, point{17}, rad2deg(motorAngle), 360*params.MeshSizeAG/params.RAG);
    mi_drawline(point{15}, point{16});
    mi_drawline(point{16}, point{17});
    mi_drawarc(point{12}, point{15}, rad2deg(motorAngle), 10);

    % BOUNDARY CONDITIONS 
    mi_addboundprop('Dirichlet',0,0,0,0);
    mi_selectarcsegment(point{12});
    mi_setarcsegmentprop(10,'Dirichlet',0,0);
    mi_clearselected

    mi_addboundprop('AP1',0,0,0,0,0,0,0,0,5);
    mi_selectsegment((point{12}+point{13})/2);
    mi_selectsegment((point{15}+point{16})/2);
    mi_setsegmentprop('AP1',0,1,0,0);
    mi_clearselected
    mi_addboundprop('AP2',0,0,0,0,0,0,0,0,5);
    mi_selectsegment((point{13}+point{14})/2);
    mi_selectsegment((point{16}+point{17})/2);
    mi_setsegmentprop('AP2',0,1,0,0);
    mi_clearselected
    
    mi_drawarc(point{13}, point{16}, rad2deg(motorAngle), 360*params.MeshSizeAG/params.RAG);
    % upper magnet
    mi_drawline(point{2}, point{6});
    mi_drawline(point{6}, point{5});
    mi_drawline(point{5}, point{1});
    mi_drawline(point{1}, point{3});
    mi_drawline(point{3}, point{9});
    mi_drawline(point{9}, point{10});
    mi_drawline(point{10}, point{4});
    mi_drawline(point{4}, point{2});
    mi_drawline(point{2}, point{6});
    mi_drawline(point{1}, point{2});
    mi_drawline(point{3}, point{4});

    % lower magnet
    mi_drawline(point{2+17}, point{6+17});
    mi_drawline(point{6+17}, point{5+17});
    mi_drawline(point{5+17}, point{1+17});
    mi_drawline(point{1+17}, point{3+17});
    mi_drawline(point{3+17}, point{9+17});
    mi_drawline(point{9+17}, point{10+17});
    mi_drawline(point{10+17}, point{4+17});
    mi_drawline(point{4+17}, point{2+17});
    mi_drawline(point{2+17}, point{6+17});
    mi_drawline(point{1+17}, point{2+17});
    mi_drawline(point{3+17}, point{4+17});
    
    % MATERIALS
    % air
    mi_getmaterial('Air');
    xm = (point{3}+point{4}+point{9}+point{10})/4;
    mi_addblocklabel(xm);
    mi_selectlabel(xm);
    xm = (point{1}+point{2}+point{5}+point{6})/4;
    mi_addblocklabel(xm);
    mi_selectlabel(xm);
    xm = (point{3+17}+point{4+17}+point{9+17}+point{10+17})/4;
    mi_addblocklabel(xm);
    mi_selectlabel(xm);
    xm = (point{1+17}+point{2+17}+point{5+17}+point{6+17})/4;
    mi_addblocklabel(xm);
    mi_selectlabel(xm);
    xm = (params.RRO+params.RAG)/2*[cos(motorAngle/2), sin(motorAngle/2)];
    mi_addblocklabel(xm);
    mi_selectlabel(xm);
    mi_setblockprop('Air', ~logical(params.MeshSize), params.MeshSize, 0, 0, 0, 0);
    mi_clearselected();

    % iron
    mi_getmaterial('M-27 Steel');
    xm = (params.RRI+params.RRO)/2*[cos(motorAngle/2), sin(motorAngle/2)];
    mi_addblocklabel(xm);
    mi_selectlabel(xm);
    mi_setblockprop('M-27 Steel', ~logical(params.MeshSize), params.MeshSize, 0, 0, 0, 0);
    mi_clearselected();

    % magnets
    mu0 = 4*pi*1e-7;
    mur = 1.05;
    Br = 1.4;
    Br = 0;
    H = Br/(mu0*mur);
    magnetization = [pi/2-deg2rad(params.MA/2)+pi/params.POLES, -pi/2+deg2rad(params.MA/2)+pi/params.POLES];
    mi_addmaterial('NdFeBr14',mur,mur,H);
    % bot magnet
    xm = (point{2}+point{3})/2;
    mi_addblocklabel(xm);
    mi_selectlabel(xm);
    mi_setblockprop('NdFeBr14', ~logical(params.MeshSize), params.MeshSize, 0, rad2deg(magnetization(2)), 0, 0);
    mi_clearselected();
    % right iron
    xm = (point{2+17}+point{3+17})/2;
    mi_addblocklabel(xm);
    mi_selectlabel(xm);
    mi_setblockprop('NdFeBr14', ~logical(params.MeshSize), params.MeshSize, 0, rad2deg(magnetization(1)), 0, 0);
    mi_clearselected();

    mi_zoomnatural();

    mi_saveas('GEO_IPMrotor_V_2.FEM');

    mytimer("FEMM-mesh", "start");
    mi_createmesh();
    mytimer("FEMM-mesh", "stop");

    mytimer("FEMM-Solve", "start");
    mi_analyze(0);
    mytimer("FEMM-Solve", "stop");

    % keyboard
    closefemm();

    mytimer("FEMM", "stop");
end