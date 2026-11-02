function params = GEO_HorseshoeMagnet_femm(params)
    arguments
        params.preview (1,1) logical = true;
        params.AirPadding = 2e-2;
        params.AirGapWidth = 2e-2;
        params.MagnetWidth = 5e-2;
        params.MagnetHeight = 2e-2;
        params.IronHeight = 5e-2; 
        params.IronWidth = 2e-2;
        params.MeshSize = 0;    % Automesh
    end 
    
    if params.MeshSize == 0
        autoMesh = true;
    else
        autoMesh = false;
    end

    widths = [0, params.AirPadding, params.IronWidth, params.MagnetWidth, params.AirGapWidth, params.IronWidth, params.AirPadding];
    cumWidths = cumsum(widths);
    heights = [0, params.AirPadding, params.MagnetHeight, params.IronHeight, params.MagnetHeight, params.AirPadding];
    cumHeights = cumsum(heights);
    
    for Hs = 1:numel(heights)
        for Ws = 1:numel(widths)
            point{Ws+(Hs-1)*numel(widths)} = [cumWidths(Ws), cumHeights(Hs)];
        end
    end

    openfemm();
    % main_minimize();
    newdocument(0);
    mi_probdef(0, 'meters', 'planar' ,1e-6, 0.1 ,30, 0);

    % GENERATE GEOMETRY
    % outer line
    mi_drawpolyline([point{1}; point{7}; point{42}; point{36}; point{1}]);
    % iron left
    mi_drawpolyline([point{9}; point{10}; point{31}; point{30}; point{9}]);
    % magnet bottom
    mi_drawpolyline([point{10}; point{11}; point{18}; point{17}]);
    % magnet top
    mi_drawpolyline([point{24}; point{25}; point{32}; point{31}]);
    % iron  right
    mi_drawpolyline([point{12}; point{13}; point{34}; point{33}; point{12}]);

    % BOUNDARY CONDITIONS 
    mi_addboundprop('Dirichlet',0,0,0,0);
    mi_selectsegment((point{1}(1)+point{7}(1))/2, (point{1}(2)+point{7}(2))/2);
    mi_selectsegment((point{7}(1)+point{42}(1))/2, (point{7}(2)+point{42}(2))/2);
    mi_selectsegment((point{42}(1)+point{36}(1))/2, (point{42}(2)+point{36}(2))/2);
    mi_selectsegment((point{36}(1)+point{1}(1))/2, (point{36}(2)+point{1}(2))/2);
    mi_setsegmentprop('Dirichlet',0,1,0,0);
    mi_clearselected

    % MATERIALS
    mi_getmaterial('Air');
    xm = (point{17}+point{25})/2;
    mi_addblocklabel(xm(1), xm(2));
    mi_selectlabel(xm(1), xm(2));
    mi_setblockprop('Air', autoMesh, params.MeshSize, 0, 0, 0, 0);
    mi_clearselected();

    mi_getmaterial('M-27 Steel');
    % left iron
    xm = (point{9}+point{31})/2;
    mi_addblocklabel(xm(1), xm(2));
    mi_selectlabel(xm(1), xm(2));
    % right iron
    xm = (point{12}+point{34})/2;
    mi_addblocklabel(xm(1), xm(2));
    mi_selectlabel(xm(1), xm(2));
    mi_setblockprop('M-27 Steel', autoMesh, params.MeshSize, 0, 0, 0, 0);
    mi_clearselected();

    mu0 = 4*pi*1e-7;
    mur = 1.05;
    Br = 1.4;
    H = Br/(mu0*mur);
    mi_addmaterial('NdFeBr14',mur,mur,H);
    % bot magnet
    xm = (point{10}+point{18})/2;
    mi_addblocklabel(xm(1), xm(2));
    mi_selectlabel(xm(1), xm(2));
    mi_setblockprop('NdFeBr14', autoMesh, params.MeshSize, 0, 180, 0, 0);
    mi_clearselected();
    % right iron
    xm = (point{24}+point{32})/2;
    mi_addblocklabel(xm(1), xm(2));
    mi_selectlabel(xm(1), xm(2));
    mi_setblockprop('NdFeBr14', autoMesh, params.MeshSize, 0, 0, 0, 0);
    mi_clearselected();

    mi_zoomnatural();

    mi_saveas('GEO_HorseshoeMagnet.FEM');

    mytimer("FEMM-mesh", "start");
    mi_createmesh();
    mytimer("FEMM-mesh", "stop");

    mytimer("FEMM-Solve", "start");
    mi_analyze(0);
    mytimer("FEMM-Solve", "stop");

    closefemm();
    
    % if (draw_geometry)
    %     figure
    %     hold on
    %     for i = 1:numel(point)
    %         scatter(point{i}(1), point{i}(2), "black", "filled");
    %         text(point{i}(1)+0.001, point{i}(2), string(i))
    %     end
    % end
end