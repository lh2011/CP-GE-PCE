% 

function params = GEO_DISstator_3_femm(params)
    arguments
        % Standard parameters   
        params.Symmetry = 8;
        params.SLOTS = 48;              % Nuber of stator slots
        params.MeshSize = 0;            % Automesh
        params.MeshSizeAG = 1e-2;        % mesh sides at air gap
    
        params.RAG = 78.9e-3;           % Radius Air Gap
        params.SRI = 79.1e-3;             % Stator Radius Inner
        params.SRO = 116e-3;            % Stator Radius Outer
        params.HSO = 1.1e-3;            % Height Slot Opening
        params.HY = 17.1e-3;            % Height Yoke
        params.WSO = 2.5e-3;            % Width Slot Opening
        params.WS = 5.1e-3;             % Width Slot
    end
    
    

    mytimer("FEMM", "start");
    
    % calculated parameters
    a1 = asin((params.WSO/2)/(params.SRI));
    a2 = pi/params.SLOTS;
    a3 = a2-a1;
    slotAngle = 2*pi/params.SLOTS;
    statorAngle = slotAngle*params.SLOTS/params.Symmetry;
    % a2 = asin((params.WSO/2)/(params.SRI+params.HSO));

    point{1} = params.RAG*[1, 0];
    point{2} = params.SRI*[1, 0];
    point{3} = params.SRO*[1, 0];
    point{4} = params.RAG*[cos(statorAngle), sin(statorAngle)];
    point{5} = params.SRI*[cos(statorAngle), sin(statorAngle)];
    point{6} = params.SRO*[cos(statorAngle), sin(statorAngle)];

    point{7} =  [cos(a1)*params.SRI, -params.WSO/2];
    LS = params.SRO-params.HY-point{7}(1);
    
    point{8} =  [cos(a1)*params.SRI+params.HSO, -params.WSO/2];
    point{9} =  [cos(a1)*params.SRI+params.HSO, -params.WS/2];
    point{10} =  point{9} + [LS, 0];
    point{11} =  point{10} + [0, params.WS];
    point{12} =  [cos(a1)*params.SRI+params.HSO, params.WS/2];
    point{13} =  [cos(a1)*params.SRI+params.HSO, params.WSO/2];
    point{14} =  [cos(a1)*params.SRI, params.WSO/2];
    point{15} = params.SRI*[cos(slotAngle), sin(slotAngle)];
    % correct unrotated points
    for iPoint = 7:14
        point{iPoint} = ([cos(slotAngle/2), -sin(slotAngle/2); sin(slotAngle/2), cos(slotAngle/2)]*point{iPoint}')';
    end
    nPoints = numel(point);
    % copy and rotate points for all slots
    for iSlot = 1:(params.SLOTS /params.Symmetry)
        for iPoint = 1:nPoints
            rotAngle = slotAngle*iSlot;
            point{nPoints+(iSlot-1)*nPoints+iPoint} =  ([cos(rotAngle), -sin(rotAngle); sin(rotAngle), cos(rotAngle)]*point{iPoint}')';
        end
    end

    % figure
    % hold on
    % for i = 1:numel(point)
    %     scatter(point{i}(1), point{i}(2), "black", "filled");
    %     text(point{i}(1)+0.0001, point{i}(2), string(i))
    % end

    openfemm();
    newdocument(0);
    mi_probdef(0, 'meters', 'planar' ,1e-6, 0.1 ,30, 0);

    % GENERATE GEOMETRY
    % outer line
    mi_drawarc(point{1}, point{4}, rad2deg(statorAngle), 360*params.MeshSizeAG/params.RAG);
    mi_drawarc(point{3}, point{6}, rad2deg(statorAngle), 1);

    mi_drawline(point{1}, point{2});
    mi_drawline(point{2}, point{3});
    mi_drawline(point{4}, point{5});
    mi_drawline(point{5}, point{6});


    % BOUNDARY CONDITIONS 
    mi_addboundprop('Dirichlet',0,0,0,0);
    mi_selectarcsegment(params.SRO*[cos(statorAngle/2), sin(statorAngle/2)]);
    mi_setarcsegmentprop(1,'Dirichlet',0,0);
    mi_clearselected

    mi_addboundprop('AP1',0,0,0,0,0,0,0,0,5);
    mi_selectsegment((point{1}+point{2})/2);
    mi_selectsegment((point{4}+point{5})/2);
    mi_setsegmentprop('AP1',0,1,0,0);
    mi_clearselected
    mi_addboundprop('AP2',0,0,0,0,0,0,0,0,5);
    mi_selectsegment((point{2}+point{3})/2);
    mi_selectsegment((point{5}+point{6})/2);
    mi_setsegmentprop('AP2',0,1,0,0);
    mi_clearselected

    
    mi_getmaterial('10 AWG');

    % create slots
    for iSlot = 1:(params.SLOTS/params.Symmetry)
        mi_drawarc(point{2+(iSlot-1)*nPoints}, point{7+(iSlot-1)*nPoints}, rad2deg(a3), 360*params.MeshSizeAG/params.SRI);
        mi_drawline(point{7+(iSlot-1)*nPoints}, point{8+(iSlot-1)*nPoints});
        mi_drawline(point{8+(iSlot-1)*nPoints}, point{9+(iSlot-1)*nPoints});
        mi_drawline(point{9+(iSlot-1)*nPoints}, point{10+(iSlot-1)*nPoints});
        mi_drawline(point{10+(iSlot-1)*nPoints}, point{11+(iSlot-1)*nPoints});
        mi_drawline(point{11+(iSlot-1)*nPoints}, point{12+(iSlot-1)*nPoints});
        mi_drawline(point{12+(iSlot-1)*nPoints}, point{13+(iSlot-1)*nPoints});
        mi_drawline(point{13+(iSlot-1)*nPoints}, point{14+(iSlot-1)*nPoints});
        mi_drawline(point{8+(iSlot-1)*nPoints}, point{13+(iSlot-1)*nPoints});
        mi_drawarc(point{14+(iSlot-1)*nPoints}, point{15+(iSlot-1)*nPoints}, rad2deg(a3), 360*params.MeshSizeAG/params.SRI);
        % set copper
        xm = (point{9+(iSlot-1)*nPoints} + point{11+(iSlot-1)*nPoints})/2;
        mi_addblocklabel(xm);
        mi_selectlabel(xm);
        mi_setblockprop('10 AWG',  ~logical(params.MeshSize), params.MeshSize, 0, 0, 0, 0);
        mi_clearselected();
    end

    mi_getmaterial('M-27 Steel');
    xm = (params.SRO-params.HY/2)*[cos(statorAngle/2), sin(statorAngle/2)];
    mi_addblocklabel(xm);
    mi_selectlabel(xm);
    mi_setblockprop('M-27 Steel',  ~logical(params.MeshSize), params.MeshSize, 0, 0, 0, 0);
    mi_clearselected();

    mi_getmaterial('Air');
    xm = (params.SRI+params.RAG)/2*[cos(statorAngle/2), sin(statorAngle/2)];
    mi_addblocklabel(xm);
    mi_selectlabel(xm);
    mi_setblockprop('Air',  ~logical(params.MeshSize), params.MeshSize, 0, 0, 0, 0);
    mi_clearselected();

    mi_zoomnatural
    
    mi_saveas('GEO_DISstator_4.FEM');

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
